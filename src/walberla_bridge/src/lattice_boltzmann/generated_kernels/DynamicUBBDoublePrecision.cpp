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
//! \\file DynamicUBBDoublePrecision.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.3.7+13.gdfd203a, lbmpy v1.3.7+5.gd7100a3, sympy v1.10, lbmpy_walberla/pystencils_walberla from waLBerla commit c69cb11d6a95d32b2280544d3d9abde1fe5fdbb5

#include "DynamicUBBDoublePrecision.h"
#include "core/DataTypes.h"
#include "core/Macros.h"

#define FUNC_PREFIX

using namespace std;

namespace walberla {
namespace lbm {

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
namespace internal_3cfabb7f34e389af9363b890d8729bed {
static FUNC_PREFIX void dynamicubbdoubleprecision_boundary_DynamicUBBDoublePrecision(uint8_t *RESTRICT _data_forceVector, uint8_t *RESTRICT const _data_indexVector, double *RESTRICT _data_pdfs, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, int32_t indexVectorSize) {

  const int32_t f_in_inv_dir_idx[] = {0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 16, 15, 18, 17, 12, 11, 14, 13};
  const int32_t f_in_inv_offsets_x[] = {0, 0, 0, -1, 1, 0, 0, -1, 1, -1, 1, 0, 0, -1, 1, 0, 0, -1, 1};
  const int32_t f_in_inv_offsets_y[] = {0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, -1, 0, 0, 1, -1, 0, 0};
  const int32_t f_in_inv_offsets_z[] = {0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1};

  const double weights[] = {((double)(0.33333333333333333)), ((double)(0.055555555555555556)), ((double)(0.055555555555555556)), ((double)(0.055555555555555556)), ((double)(0.055555555555555556)), ((double)(0.055555555555555556)), ((double)(0.055555555555555556)), ((double)(0.027777777777777778)), ((double)(0.027777777777777778)), ((double)(0.027777777777777778)), ((double)(0.027777777777777778)), ((double)(0.027777777777777778)), ((double)(0.027777777777777778)), ((double)(0.027777777777777778)), ((double)(0.027777777777777778)), ((double)(0.027777777777777778)), ((double)(0.027777777777777778)), ((double)(0.027777777777777778)), ((double)(0.027777777777777778))};

  const int32_t neighbour_offset_x[] = {0, 0, 0, -1, 1, 0, 0, -1, 1, -1, 1, 0, 0, -1, 1, 0, 0, -1, 1};
  const int32_t neighbour_offset_y[] = {0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, -1, 0, 0, 1, -1, 0, 0};
  const int32_t neighbour_offset_z[] = {0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1};

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int64_t ctr_0 = 0; ctr_0 < indexVectorSize; ctr_0 += 1) {
      const int32_t x = *((int32_t *)(&_data_indexVector[40 * ctr_0]));
      const int32_t y = *((int32_t *)(&_data_indexVector[40 * ctr_0 + 4]));
      const int32_t z = *((int32_t *)(&_data_indexVector[40 * ctr_0 + 8]));
      const int32_t dir = *((int32_t *)(&_data_indexVector[40 * ctr_0 + 12]));
      const double vel0Term = _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 10 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 14 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 18 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 4 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 8 * _stride_pdfs_3];
      const double vel1Term = _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 11 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 15 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 7 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y + _stride_pdfs_2 * z + _stride_pdfs_3];
      const double vel2Term = _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 12 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 13 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 5 * _stride_pdfs_3];
      const double delta_rho = vel0Term + vel1Term + vel2Term + _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 16 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 17 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 2 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 3 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 6 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y + _stride_pdfs_2 * z + 9 * _stride_pdfs_3] + _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y + _stride_pdfs_2 * z];
      const double rho = delta_rho + 1.0;
      _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_0 * f_in_inv_offsets_x[dir] + _stride_pdfs_1 * y + _stride_pdfs_1 * f_in_inv_offsets_y[dir] + _stride_pdfs_2 * z + _stride_pdfs_2 * f_in_inv_offsets_z[dir] + _stride_pdfs_3 * f_in_inv_dir_idx[dir]] = -rho * (6.0 * ((double)(neighbour_offset_x[dir])) * *((double *)(&_data_indexVector[40 * ctr_0 + 16])) + 6.0 * ((double)(neighbour_offset_y[dir])) * *((double *)(&_data_indexVector[40 * ctr_0 + 24])) + 6.0 * ((double)(neighbour_offset_z[dir])) * *((double *)(&_data_indexVector[40 * ctr_0 + 32]))) * weights[dir] + _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y + _stride_pdfs_2 * z + _stride_pdfs_3 * dir];
      const double f = -rho * (6.0 * ((double)(neighbour_offset_x[dir])) * *((double *)(&_data_indexVector[40 * ctr_0 + 16])) + 6.0 * ((double)(neighbour_offset_y[dir])) * *((double *)(&_data_indexVector[40 * ctr_0 + 24])) + 6.0 * ((double)(neighbour_offset_z[dir])) * *((double *)(&_data_indexVector[40 * ctr_0 + 32]))) * weights[dir] + 2.0 * _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y + _stride_pdfs_2 * z + _stride_pdfs_3 * dir];
      *((double *)(&_data_forceVector[24 * ctr_0])) = f * ((double)(neighbour_offset_x[dir]));
      *((double *)(&_data_forceVector[24 * ctr_0 + 8])) = f * ((double)(neighbour_offset_y[dir]));
      *((double *)(&_data_forceVector[24 * ctr_0 + 16])) = f * ((double)(neighbour_offset_z[dir]));
    }
  }
}
} // namespace internal_3cfabb7f34e389af9363b890d8729bed

// NOLINTEND(readability-non-const-parameter*)
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#ifdef __CUDACC__
#pragma pop
#endif

void DynamicUBBDoublePrecision::run_impl(IBlock *block, IndexVectors::Type type) {
  auto *indexVectors = block->getData<IndexVectors>(indexVectorID);
  int32_t indexVectorSize = int32_c(indexVectors->indexVector(type).size());
  if (indexVectorSize == 0)
    return;

  auto pointer = indexVectors->pointerCpu(type);

  uint8_t *_data_indexVector = reinterpret_cast<uint8_t *>(pointer);

  auto *forceVector = block->getData<ForceVector>(forceVectorID);
  WALBERLA_ASSERT_EQUAL(indexVectorSize, int32_c(forceVector->forceVector().size()))

  auto forcePointer = forceVector->pointerCpu();

  uint8_t *_data_forceVector = reinterpret_cast<uint8_t *>(forcePointer);

  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);

  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()))
  double *RESTRICT _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
  const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  internal_3cfabb7f34e389af9363b890d8729bed::dynamicubbdoubleprecision_boundary_DynamicUBBDoublePrecision(_data_forceVector, _data_indexVector, _data_pdfs, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, indexVectorSize);
}

void DynamicUBBDoublePrecision::run(IBlock *block) {
  run_impl(block, IndexVectors::ALL);
}

void DynamicUBBDoublePrecision::inner(IBlock *block) {
  run_impl(block, IndexVectors::INNER);
}

void DynamicUBBDoublePrecision::outer(IBlock *block) {
  run_impl(block, IndexVectors::OUTER);
}

} // namespace lbm
} // namespace walberla
