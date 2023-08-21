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
//! \\file ReactionKernelBulk_4_double_precision.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

// kernel generated with pystencils v1.2, lbmpy v1.3.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 065ce5f311850371a97ac4766f47dbb5ca8424ba

#include <cmath>

#include "ReactionKernelBulk_4_double_precision.h"
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

namespace internal_49f1aaa6789b7fa16fb103a21ce6fe12 {
static FUNC_PREFIX void reactionkernelbulk_4_double_precision_reactionkernelbulk_4_double_precision(double *RESTRICT _data_rho_0, double *RESTRICT _data_rho_1, double *RESTRICT _data_rho_2, double *RESTRICT _data_rho_3, int64_t const _size_rho_0_0, int64_t const _size_rho_0_1, int64_t const _size_rho_0_2, int64_t const _stride_rho_0_0, int64_t const _stride_rho_0_1, int64_t const _stride_rho_0_2, int64_t const _stride_rho_1_0, int64_t const _stride_rho_1_1, int64_t const _stride_rho_1_2, int64_t const _stride_rho_2_0, int64_t const _stride_rho_2_1, int64_t const _stride_rho_2_2, int64_t const _stride_rho_3_0, int64_t const _stride_rho_3_1, int64_t const _stride_rho_3_2, double order_0, double order_1, double order_2, double order_3, double rate_coefficient, double stoech_0, double stoech_1, double stoech_2, double stoech_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_rho_0_2; ctr_2 += 1) {
    double *RESTRICT _data_rho_0_20 = _data_rho_0 + _stride_rho_0_2 * ctr_2;
    double *RESTRICT _data_rho_1_20 = _data_rho_1 + _stride_rho_1_2 * ctr_2;
    double *RESTRICT _data_rho_2_20 = _data_rho_2 + _stride_rho_2_2 * ctr_2;
    double *RESTRICT _data_rho_3_20 = _data_rho_3 + _stride_rho_3_2 * ctr_2;
    for (int64_t ctr_1 = 0; ctr_1 < _size_rho_0_1; ctr_1 += 1) {
      double *RESTRICT _data_rho_0_20_10 = _stride_rho_0_1 * ctr_1 + _data_rho_0_20;
      double *RESTRICT _data_rho_1_20_10 = _stride_rho_1_1 * ctr_1 + _data_rho_1_20;
      double *RESTRICT _data_rho_2_20_10 = _stride_rho_2_1 * ctr_1 + _data_rho_2_20;
      double *RESTRICT _data_rho_3_20_10 = _stride_rho_3_1 * ctr_1 + _data_rho_3_20;
      for (int64_t ctr_0 = 0; ctr_0 < _size_rho_0_0; ctr_0 += 1) {
        const double local_rho_0 = _data_rho_0_20_10[_stride_rho_0_0 * ctr_0];
        const double local_rho_1 = _data_rho_1_20_10[_stride_rho_1_0 * ctr_0];
        const double local_rho_2 = _data_rho_2_20_10[_stride_rho_2_0 * ctr_0];
        const double local_rho_3 = _data_rho_3_20_10[_stride_rho_3_0 * ctr_0];
        const double rate_factor = pow(local_rho_0, order_0) * pow(local_rho_1, order_1) * pow(local_rho_2, order_2) * pow(local_rho_3, order_3) * rate_coefficient;
        _data_rho_0_20_10[_stride_rho_0_0 * ctr_0] = local_rho_0 + rate_factor * stoech_0;
        _data_rho_1_20_10[_stride_rho_1_0 * ctr_0] = local_rho_1 + rate_factor * stoech_1;
        _data_rho_2_20_10[_stride_rho_2_0 * ctr_0] = local_rho_2 + rate_factor * stoech_2;
        _data_rho_3_20_10[_stride_rho_3_0 * ctr_0] = local_rho_3 + rate_factor * stoech_3;
      }
    }
  }
}
} // namespace internal_49f1aaa6789b7fa16fb103a21ce6fe12

void ReactionKernelBulk_4_double_precision::run(IBlock *block) {
  auto rho_0 = block->getData<field::GhostLayerField<double, 1>>(rho_0ID);
  auto rho_2 = block->getData<field::GhostLayerField<double, 1>>(rho_2ID);
  auto rho_1 = block->getData<field::GhostLayerField<double, 1>>(rho_1ID);
  auto rho_3 = block->getData<field::GhostLayerField<double, 1>>(rho_3ID);

  auto &order_1 = this->order_1_;
  auto &order_2 = this->order_2_;
  auto &stoech_0 = this->stoech_0_;
  auto &stoech_2 = this->stoech_2_;
  auto &order_3 = this->order_3_;
  auto &stoech_1 = this->stoech_1_;
  auto &stoech_3 = this->stoech_3_;
  auto &order_0 = this->order_0_;
  auto &rate_coefficient = this->rate_coefficient_;
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(rho_0->nrOfGhostLayers()));
  double *RESTRICT _data_rho_0 = rho_0->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(rho_1->nrOfGhostLayers()));
  double *RESTRICT _data_rho_1 = rho_1->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(rho_2->nrOfGhostLayers()));
  double *RESTRICT _data_rho_2 = rho_2->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(rho_3->nrOfGhostLayers()));
  double *RESTRICT _data_rho_3 = rho_3->dataAt(0, 0, 0, 0);
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
  internal_49f1aaa6789b7fa16fb103a21ce6fe12::reactionkernelbulk_4_double_precision_reactionkernelbulk_4_double_precision(_data_rho_0, _data_rho_1, _data_rho_2, _data_rho_3, _size_rho_0_0, _size_rho_0_1, _size_rho_0_2, _stride_rho_0_0, _stride_rho_0_1, _stride_rho_0_2, _stride_rho_1_0, _stride_rho_1_1, _stride_rho_1_2, _stride_rho_2_0, _stride_rho_2_1, _stride_rho_2_2, _stride_rho_3_0, _stride_rho_3_1, _stride_rho_3_2, order_0, order_1, order_2, order_3, rate_coefficient, stoech_0, stoech_1, stoech_2, stoech_3);
}

void ReactionKernelBulk_4_double_precision::runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks, const CellInterval &globalCellInterval, cell_idx_t ghostLayers, IBlock *block) {
  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto rho_0 = block->getData<field::GhostLayerField<double, 1>>(rho_0ID);
  auto rho_2 = block->getData<field::GhostLayerField<double, 1>>(rho_2ID);
  auto rho_1 = block->getData<field::GhostLayerField<double, 1>>(rho_1ID);
  auto rho_3 = block->getData<field::GhostLayerField<double, 1>>(rho_3ID);

  auto &order_1 = this->order_1_;
  auto &order_2 = this->order_2_;
  auto &stoech_0 = this->stoech_0_;
  auto &stoech_2 = this->stoech_2_;
  auto &order_3 = this->order_3_;
  auto &stoech_1 = this->stoech_1_;
  auto &stoech_3 = this->stoech_3_;
  auto &order_0 = this->order_0_;
  auto &rate_coefficient = this->rate_coefficient_;
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(rho_0->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(rho_0->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(rho_0->nrOfGhostLayers()));
  double *RESTRICT _data_rho_0 = rho_0->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(rho_1->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(rho_1->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(rho_1->nrOfGhostLayers()));
  double *RESTRICT _data_rho_1 = rho_1->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(rho_2->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(rho_2->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(rho_2->nrOfGhostLayers()));
  double *RESTRICT _data_rho_2 = rho_2->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(rho_3->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(rho_3->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(rho_3->nrOfGhostLayers()));
  double *RESTRICT _data_rho_3 = rho_3->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
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
  internal_49f1aaa6789b7fa16fb103a21ce6fe12::reactionkernelbulk_4_double_precision_reactionkernelbulk_4_double_precision(_data_rho_0, _data_rho_1, _data_rho_2, _data_rho_3, _size_rho_0_0, _size_rho_0_1, _size_rho_0_2, _stride_rho_0_0, _stride_rho_0_1, _stride_rho_0_2, _stride_rho_1_0, _stride_rho_1_1, _stride_rho_1_2, _stride_rho_2_0, _stride_rho_2_1, _stride_rho_2_2, _stride_rho_3_0, _stride_rho_3_1, _stride_rho_3_2, order_0, order_1, order_2, order_3, rate_coefficient, stoech_0, stoech_1, stoech_2, stoech_3);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif