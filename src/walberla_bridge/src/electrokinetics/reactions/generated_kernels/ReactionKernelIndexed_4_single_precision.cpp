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
//! \\file ReactionKernelIndexed_4_single_precision.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.3.3, lbmpy v1.3.3, lbmpy_walberla/pystencils_walberla from waLBerla commit b0842e1a493ce19ef1bbb8d2cf382fc343970a7f

#include "ReactionKernelIndexed_4_single_precision.h"
#include "core/DataTypes.h"
#include "core/Macros.h"

#define FUNC_PREFIX

using namespace std;

namespace walberla {
namespace pystencils {

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wconversion"
#endif

#ifdef __CUDACC__
#pragma push
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diag_suppress 177
#else
#pragma diag_suppress 177
#endif
#endif
// NOLINTBEGIN(readability-non-const-parameter*)
namespace internal_62195c7636fc531694c87493561447d7 {
static FUNC_PREFIX void reactionkernelindexed_4_single_precision_boundary_ReactionKernelIndexed_4_single_precision(uint8_t *RESTRICT const _data_indexVector, float *RESTRICT _data_rho_0, float *RESTRICT _data_rho_1, float *RESTRICT _data_rho_2, float *RESTRICT _data_rho_3, int64_t const _stride_rho_0_0, int64_t const _stride_rho_0_1, int64_t const _stride_rho_0_2, int64_t const _stride_rho_1_0, int64_t const _stride_rho_1_1, int64_t const _stride_rho_1_2, int64_t const _stride_rho_2_0, int64_t const _stride_rho_2_1, int64_t const _stride_rho_2_2, int64_t const _stride_rho_3_0, int64_t const _stride_rho_3_1, int64_t const _stride_rho_3_2, int32_t indexVectorSize, float order_0, float order_1, float order_2, float order_3, float rate_coefficient, float stoech_0, float stoech_1, float stoech_2, float stoech_3) {
  for (int64_t ctr_0 = 0; ctr_0 < indexVectorSize; ctr_0 += 1) {
    const int32_t x = *((int32_t *)(&_data_indexVector[12 * ctr_0]));
    const int32_t y = *((int32_t *)(&_data_indexVector[12 * ctr_0 + 4]));
    const int32_t z = *((int32_t *)(&_data_indexVector[12 * ctr_0 + 8]));

    const int32_t cx[] = {0};
    const int32_t cy[] = {0};
    const int32_t cz[] = {0};
    const int32_t invdir[] = {0};

    const float local_rho_0 = _data_rho_0[_stride_rho_0_0 * x + _stride_rho_0_1 * y + _stride_rho_0_2 * z];
    const float local_rho_1 = _data_rho_1[_stride_rho_1_0 * x + _stride_rho_1_1 * y + _stride_rho_1_2 * z];
    const float local_rho_2 = _data_rho_2[_stride_rho_2_0 * x + _stride_rho_2_1 * y + _stride_rho_2_2 * z];
    const float local_rho_3 = _data_rho_3[_stride_rho_3_0 * x + _stride_rho_3_1 * y + _stride_rho_3_2 * z];
    const float rate_factor = rate_coefficient * powf(local_rho_0, order_0) * powf(local_rho_1, order_1) * powf(local_rho_2, order_2) * powf(local_rho_3, order_3);
    _data_rho_0[_stride_rho_0_0 * x + _stride_rho_0_1 * y + _stride_rho_0_2 * z] = local_rho_0 + rate_factor * stoech_0;
    _data_rho_1[_stride_rho_1_0 * x + _stride_rho_1_1 * y + _stride_rho_1_2 * z] = local_rho_1 + rate_factor * stoech_1;
    _data_rho_2[_stride_rho_2_0 * x + _stride_rho_2_1 * y + _stride_rho_2_2 * z] = local_rho_2 + rate_factor * stoech_2;
    _data_rho_3[_stride_rho_3_0 * x + _stride_rho_3_1 * y + _stride_rho_3_2 * z] = local_rho_3 + rate_factor * stoech_3;
  }
}
} // namespace internal_62195c7636fc531694c87493561447d7

// NOLINTEND(readability-non-const-parameter*)
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#ifdef __CUDACC__
#pragma pop
#endif

void ReactionKernelIndexed_4_single_precision::run_impl(IBlock *block, IndexVectors::Type type) {
  auto *indexVectors = block->uncheckedFastGetData<IndexVectors>(indexVectorID);
  int32_t indexVectorSize = int32_c(indexVectors->indexVector(type).size());
  if (indexVectorSize == 0)
    return;

  auto pointer = indexVectors->pointerCpu(type);

  uint8_t *_data_indexVector = reinterpret_cast<uint8_t *>(pointer);

  auto rho_3 = block->getData<field::GhostLayerField<float, 1>>(rho_3ID);
  auto rho_1 = block->getData<field::GhostLayerField<float, 1>>(rho_1ID);
  auto rho_2 = block->getData<field::GhostLayerField<float, 1>>(rho_2ID);
  auto rho_0 = block->getData<field::GhostLayerField<float, 1>>(rho_0ID);

  auto &order_1 = order_1_;
  auto &order_0 = order_0_;
  auto &order_2 = order_2_;
  auto &stoech_0 = stoech_0_;
  auto &stoech_1 = stoech_1_;
  auto &order_3 = order_3_;
  auto &stoech_2 = stoech_2_;
  auto &rate_coefficient = rate_coefficient_;
  auto &stoech_3 = stoech_3_;
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(rho_0->nrOfGhostLayers()))
  float *RESTRICT _data_rho_0 = rho_0->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(rho_1->nrOfGhostLayers()))
  float *RESTRICT _data_rho_1 = rho_1->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(rho_2->nrOfGhostLayers()))
  float *RESTRICT _data_rho_2 = rho_2->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(rho_3->nrOfGhostLayers()))
  float *RESTRICT _data_rho_3 = rho_3->dataAt(0, 0, 0, 0);
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
  internal_62195c7636fc531694c87493561447d7::reactionkernelindexed_4_single_precision_boundary_ReactionKernelIndexed_4_single_precision(_data_indexVector, _data_rho_0, _data_rho_1, _data_rho_2, _data_rho_3, _stride_rho_0_0, _stride_rho_0_1, _stride_rho_0_2, _stride_rho_1_0, _stride_rho_1_1, _stride_rho_1_2, _stride_rho_2_0, _stride_rho_2_1, _stride_rho_2_2, _stride_rho_3_0, _stride_rho_3_1, _stride_rho_3_2, indexVectorSize, order_0, order_1, order_2, order_3, rate_coefficient, stoech_0, stoech_1, stoech_2, stoech_3);
}

void ReactionKernelIndexed_4_single_precision::run(IBlock *block) {
  run_impl(block, IndexVectors::ALL);
}

void ReactionKernelIndexed_4_single_precision::inner(IBlock *block) {
  run_impl(block, IndexVectors::INNER);
}

void ReactionKernelIndexed_4_single_precision::outer(IBlock *block) {
  run_impl(block, IndexVectors::OUTER);
}

} // namespace pystencils
} // namespace walberla
