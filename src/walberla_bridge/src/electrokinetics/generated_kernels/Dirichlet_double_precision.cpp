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
//! \\file Dirichlet_double_precision.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.2, lbmpy v1.2, lbmpy_walberla/pystencils_walberla from waLBerla commit ref: a839fac6ef7d0c58e7710e4d50490e9dd7146b4a

#include <cmath>

#include "Dirichlet_double_precision.h"
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

namespace internal_74da74b67a122b7887d3d21c7ea5f414 {
static FUNC_PREFIX void dirichlet_double_precision_boundary_Dirichlet_double_precision(double *RESTRICT _data_field, uint8_t *RESTRICT const _data_indexVector, int64_t const _stride_field_0, int64_t const _stride_field_1, int64_t const _stride_field_2, int32_t indexVectorSize) {
  for (int64_t ctr_0 = 0; ctr_0 < indexVectorSize; ctr_0 += 1) {
    const int32_t x = *((int32_t *)(&_data_indexVector[24 * ctr_0]));
    const int32_t y = *((int32_t *)(&_data_indexVector[24 * ctr_0 + 4]));
    const int32_t z = *((int32_t *)(&_data_indexVector[24 * ctr_0 + 8]));

    const int32_t cx[] = {0, 0, 0, -1, 1, 0, 0, -1, 1, -1, 1, 0, 0, -1, 1, 0, 0, -1, 1, 1, -1, 1, -1, 1, -1, 1, -1};
    const int32_t cy[] = {0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, -1, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1};
    const int32_t cz[] = {0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1};
    const int32_t invdir[] = {0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 16, 15, 18, 17, 12, 11, 14, 13, 26, 25, 24, 23, 22, 21, 20, 19};

    const int32_t dir = *((int32_t *)(&_data_indexVector[24 * ctr_0 + 12]));
    _data_field[_stride_field_0 * x + _stride_field_1 * y + _stride_field_2 * z] = *((double *)(&_data_indexVector[24 * ctr_0 + 16]));
  }
}
} // namespace internal_74da74b67a122b7887d3d21c7ea5f414

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#ifdef __CUDACC__
#pragma pop
#endif

void Dirichlet_double_precision::run_impl(IBlock *block, IndexVectors::Type type) {
  auto *indexVectors = block->getData<IndexVectors>(indexVectorID);
  int32_t indexVectorSize = int32_c(indexVectors->indexVector(type).size());
  if (indexVectorSize == 0)
    return;

  auto pointer = indexVectors->pointerCpu(type);

  uint8_t *_data_indexVector = reinterpret_cast<uint8_t *>(pointer);

  auto field = block->getData<field::GhostLayerField<double, 1>>(fieldID);

  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(field->nrOfGhostLayers()));
  double *RESTRICT _data_field = field->dataAt(0, 0, 0, 0);
  const int64_t _stride_field_0 = int64_t(field->xStride());
  const int64_t _stride_field_1 = int64_t(field->yStride());
  const int64_t _stride_field_2 = int64_t(field->zStride());
  internal_74da74b67a122b7887d3d21c7ea5f414::dirichlet_double_precision_boundary_Dirichlet_double_precision(_data_field, _data_indexVector, _stride_field_0, _stride_field_1, _stride_field_2, indexVectorSize);
}

void Dirichlet_double_precision::run(IBlock *block) {
  run_impl(block, IndexVectors::ALL);
}

void Dirichlet_double_precision::inner(IBlock *block) {
  run_impl(block, IndexVectors::INNER);
}

void Dirichlet_double_precision::outer(IBlock *block) {
  run_impl(block, IndexVectors::OUTER);
}

} // namespace pystencils
} // namespace walberla
