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
//! \\file ReactionKernelBulk_1_double_precision.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.3.3, lbmpy v1.3.3, lbmpy_walberla/pystencils_walberla from waLBerla commit b0842e1a493ce19ef1bbb8d2cf382fc343970a7f

#include <cmath>

#include "ReactionKernelBulk_1_double_precision.h"
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

namespace internal_a94c3c474646ee6905a4b90e8ccc47e6 {
static FUNC_PREFIX void reactionkernelbulk_1_double_precision_reactionkernelbulk_1_double_precision(double *RESTRICT _data_rho_0, int64_t const _size_rho_0_0, int64_t const _size_rho_0_1, int64_t const _size_rho_0_2, int64_t const _stride_rho_0_0, int64_t const _stride_rho_0_1, int64_t const _stride_rho_0_2, double order_0, double rate_coefficient, double stoech_0) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_rho_0_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_rho_0_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_rho_0_0; ctr_0 += 1) {
        const double local_rho_0 = _data_rho_0[_stride_rho_0_0 * ctr_0 + _stride_rho_0_1 * ctr_1 + _stride_rho_0_2 * ctr_2];
        const double rate_factor = pow(local_rho_0, order_0) * rate_coefficient;
        _data_rho_0[_stride_rho_0_0 * ctr_0 + _stride_rho_0_1 * ctr_1 + _stride_rho_0_2 * ctr_2] = local_rho_0 + rate_factor * stoech_0;
      }
    }
  }
}
} // namespace internal_a94c3c474646ee6905a4b90e8ccc47e6

void ReactionKernelBulk_1_double_precision::run(IBlock *block) {

  auto rho_0 = block->getData<field::GhostLayerField<double, 1>>(rho_0ID);

  auto &stoech_0 = this->stoech_0_;
  auto &rate_coefficient = this->rate_coefficient_;
  auto &order_0 = this->order_0_;
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(rho_0->nrOfGhostLayers()))
  double *RESTRICT _data_rho_0 = rho_0->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(rho_0->xSizeWithGhostLayer(), int64_t(int64_c(rho_0->xSize()) + 0))
  const int64_t _size_rho_0_0 = int64_t(int64_c(rho_0->xSize()) + 0);
  WALBERLA_ASSERT_GREATER_EQUAL(rho_0->ySizeWithGhostLayer(), int64_t(int64_c(rho_0->ySize()) + 0))
  const int64_t _size_rho_0_1 = int64_t(int64_c(rho_0->ySize()) + 0);
  WALBERLA_ASSERT_GREATER_EQUAL(rho_0->zSizeWithGhostLayer(), int64_t(int64_c(rho_0->zSize()) + 0))
  const int64_t _size_rho_0_2 = int64_t(int64_c(rho_0->zSize()) + 0);
  const int64_t _stride_rho_0_0 = int64_t(rho_0->xStride());
  const int64_t _stride_rho_0_1 = int64_t(rho_0->yStride());
  const int64_t _stride_rho_0_2 = int64_t(rho_0->zStride());
  internal_a94c3c474646ee6905a4b90e8ccc47e6::reactionkernelbulk_1_double_precision_reactionkernelbulk_1_double_precision(_data_rho_0, _size_rho_0_0, _size_rho_0_1, _size_rho_0_2, _stride_rho_0_0, _stride_rho_0_1, _stride_rho_0_2, order_0, rate_coefficient, stoech_0);
}

void ReactionKernelBulk_1_double_precision::runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks, const CellInterval &globalCellInterval, cell_idx_t ghostLayers, IBlock *block) {

  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto rho_0 = block->getData<field::GhostLayerField<double, 1>>(rho_0ID);

  auto &stoech_0 = this->stoech_0_;
  auto &rate_coefficient = this->rate_coefficient_;
  auto &order_0 = this->order_0_;
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(rho_0->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(rho_0->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(rho_0->nrOfGhostLayers()))
  double *RESTRICT _data_rho_0 = rho_0->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_GREATER_EQUAL(rho_0->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
  const int64_t _size_rho_0_0 = int64_t(int64_c(ci.xSize()) + 0);
  WALBERLA_ASSERT_GREATER_EQUAL(rho_0->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
  const int64_t _size_rho_0_1 = int64_t(int64_c(ci.ySize()) + 0);
  WALBERLA_ASSERT_GREATER_EQUAL(rho_0->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
  const int64_t _size_rho_0_2 = int64_t(int64_c(ci.zSize()) + 0);
  const int64_t _stride_rho_0_0 = int64_t(rho_0->xStride());
  const int64_t _stride_rho_0_1 = int64_t(rho_0->yStride());
  const int64_t _stride_rho_0_2 = int64_t(rho_0->zStride());
  internal_a94c3c474646ee6905a4b90e8ccc47e6::reactionkernelbulk_1_double_precision_reactionkernelbulk_1_double_precision(_data_rho_0, _size_rho_0_0, _size_rho_0_1, _size_rho_0_2, _stride_rho_0_0, _stride_rho_0_1, _stride_rho_0_2, order_0, rate_coefficient, stoech_0);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif
