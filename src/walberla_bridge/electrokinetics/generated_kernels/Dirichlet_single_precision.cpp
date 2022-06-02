// kernel generated with pystencils v0.4.4+2.g825be1d, lbmpy v0.4.4,
// lbmpy_walberla/pystencils_walberla from commit
// 3c4ea5df560c4dfdd15076648fc59c830bd64a9d

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
//! \\file Dirichlet_single_precision.cpp
//! \\author pystencils
//======================================================================================================================

#include <cmath>

#include "Dirichlet_single_precision.h"
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

namespace internal_43f4eae176e72ad2d9db0f0468064c30 {
static FUNC_PREFIX void
dirichlet_single_precision_boundary_Dirichlet_single_precision(
    float *RESTRICT _data_field, uint8_t *RESTRICT const _data_indexVector,
    int64_t const _stride_field_0, int64_t const _stride_field_1,
    int64_t const _stride_field_2, int64_t indexVectorSize) {
  for (int64_t ctr_0 = 0; ctr_0 < indexVectorSize; ctr_0 += 1) {
    const int64_t x = *((int32_t *)(&_data_indexVector[20 * ctr_0]));
    const int64_t y = *((int32_t *)(&_data_indexVector[20 * ctr_0 + 4]));
    const int64_t z = *((int32_t *)(&_data_indexVector[20 * ctr_0 + 8]));

    const int64_t cx[] = {0, 0, 0, -1, 1, 0, 0,  -1, 1,  -1, 1,  0, 0, -1,
                          1, 0, 0, -1, 1, 1, -1, 1,  -1, 1,  -1, 1, -1};
    const int64_t cy[] = {0, 1, -1, 0, 0, 0, 0, 1,  1,  -1, -1, 1,  -1, 0,
                          0, 1, -1, 0, 0, 1, 1, -1, -1, 1,  1,  -1, -1};
    const int64_t cz[] = {0, 0,  0,  0,  0,  1, -1, 0, 0, 0,  0,  1,  1, 1,
                          1, -1, -1, -1, -1, 1, 1,  1, 1, -1, -1, -1, -1};
    const int64_t invdir[] = {0,  2,  1,  4,  3,  6,  5,  10, 9,
                              8,  7,  16, 15, 18, 17, 12, 11, 14,
                              13, 26, 25, 24, 23, 22, 21, 20, 19};

    const int64_t dir = *((int32_t *)(&_data_indexVector[20 * ctr_0 + 12]));
    _data_field[_stride_field_0 * x + _stride_field_1 * y +
                _stride_field_2 * z] =
        *((float *)(&_data_indexVector[20 * ctr_0 + 16]));
  }
}
} // namespace internal_43f4eae176e72ad2d9db0f0468064c30

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#ifdef __CUDACC__
#pragma pop
#endif

void Dirichlet_single_precision::run_impl(IBlock *block,
                                          IndexVectors::Type type) {
  auto *indexVectors = block->uncheckedFastGetData<IndexVectors>(indexVectorID);
  int64_t indexVectorSize = int64_c(indexVectors->indexVector(type).size());
  if (indexVectorSize == 0)
    return;

  auto pointer = indexVectors->pointerCpu(type);

  uint8_t *_data_indexVector = reinterpret_cast<uint8_t *>(pointer);

  auto field = block->getData<field::GhostLayerField<float, 1>>(fieldID);

  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(field->nrOfGhostLayers()));
  float *RESTRICT _data_field = field->dataAt(0, 0, 0, 0);
  const int64_t _stride_field_0 = int64_t(field->xStride());
  const int64_t _stride_field_1 = int64_t(field->yStride());
  const int64_t _stride_field_2 = int64_t(field->zStride());
  internal_43f4eae176e72ad2d9db0f0468064c30::
      dirichlet_single_precision_boundary_Dirichlet_single_precision(
          _data_field, _data_indexVector, _stride_field_0, _stride_field_1,
          _stride_field_2, indexVectorSize);
}

void Dirichlet_single_precision::run(IBlock *block) {
  run_impl(block, IndexVectors::ALL);
}

void Dirichlet_single_precision::inner(IBlock *block) {
  run_impl(block, IndexVectors::INNER);
}

void Dirichlet_single_precision::outer(IBlock *block) {
  run_impl(block, IndexVectors::OUTER);
}

} // namespace pystencils
} // namespace walberla
