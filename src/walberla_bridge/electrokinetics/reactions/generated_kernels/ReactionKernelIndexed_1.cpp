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
//! \\file ReactionKernelIndexed_1.cpp
//! \\author pystencils
//======================================================================================================================

#include <cmath>

#include "ReactionKernelIndexed_1.h"
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
#pragma diag_suppress = declared_but_not_referenced
#endif

namespace internal_reactionkernelindexed_1_boundary_ReactionKernelIndexed_1 {
static FUNC_PREFIX void
reactionkernelindexed_1_boundary_ReactionKernelIndexed_1(
    uint8_t *RESTRICT _data_indexVector, double *RESTRICT _data_rho_0,
    int64_t const _stride_rho_0_0, int64_t const _stride_rho_0_1,
    int64_t const _stride_rho_0_2, int64_t indexVectorSize, double order_0,
    double rate_coefficient, double stoech_0) {
  for (int64_t ctr_0 = 0; ctr_0 < indexVectorSize; ctr_0 += 1) {
    const int64_t x = *((int32_t *)(&_data_indexVector[12 * ctr_0]));
    const int64_t y = *((int32_t *)(&_data_indexVector[12 * ctr_0 + 4]));
    const int64_t z = *((int32_t *)(&_data_indexVector[12 * ctr_0 + 8]));
    const double local_rho_0 =
        _data_rho_0[_stride_rho_0_0 * x + _stride_rho_0_1 * y +
                    _stride_rho_0_2 * z];
    const double rate_factor = pow(local_rho_0, order_0) * rate_coefficient;
    _data_rho_0[_stride_rho_0_0 * x + _stride_rho_0_1 * y +
                _stride_rho_0_2 * z] = local_rho_0 + rate_factor * stoech_0;
  }
}
} // namespace internal_reactionkernelindexed_1_boundary_ReactionKernelIndexed_1

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#ifdef __CUDACC__
#pragma pop
#endif

void ReactionKernelIndexed_1::run_impl(IBlock *block, IndexVectors::Type type) {
  auto *indexVectors = block->uncheckedFastGetData<IndexVectors>(indexVectorID);
  int64_t indexVectorSize = int64_c(indexVectors->indexVector(type).size());
  if (indexVectorSize == 0)
    return;

  auto pointer = indexVectors->pointerCpu(type);

  uint8_t *_data_indexVector = reinterpret_cast<uint8_t *>(pointer);

  auto rho_0 = block->getData<field::GhostLayerField<double, 1>>(rho_0ID);

  auto &order_0 = order_0_;
  auto &rate_coefficient = rate_coefficient_;
  auto &stoech_0 = stoech_0_;
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(rho_0->nrOfGhostLayers()));
  double *RESTRICT _data_rho_0 = rho_0->dataAt(0, 0, 0, 0);
  const int64_t _stride_rho_0_0 = int64_t(rho_0->xStride());
  const int64_t _stride_rho_0_1 = int64_t(rho_0->yStride());
  const int64_t _stride_rho_0_2 = int64_t(rho_0->zStride());
  internal_reactionkernelindexed_1_boundary_ReactionKernelIndexed_1::
      reactionkernelindexed_1_boundary_ReactionKernelIndexed_1(
          _data_indexVector, _data_rho_0, _stride_rho_0_0, _stride_rho_0_1,
          _stride_rho_0_2, indexVectorSize, order_0, rate_coefficient,
          stoech_0);
}

void ReactionKernelIndexed_1::run(IBlock *block) {
  run_impl(block, IndexVectors::ALL);
}

void ReactionKernelIndexed_1::inner(IBlock *block) {
  run_impl(block, IndexVectors::INNER);
}

void ReactionKernelIndexed_1::outer(IBlock *block) {
  run_impl(block, IndexVectors::OUTER);
}

} // namespace pystencils
} // namespace walberla
