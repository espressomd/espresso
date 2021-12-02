// kernel generated with pystencils v0.4.3+12.g29e0e84, lbmpy v0.4.3+2.g0e17e61,
// lbmpy_walberla/pystencils_walberla from commit
// 08f04ef64f95609b47838db85862033a1600afa1

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
//! \\file Dynamic_UBB_double_precision.cpp
//! \\author pystencils
//======================================================================================================================

#include <cmath>

#include "Dynamic_UBB_double_precision.h"
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
#pragma diag_suppress = declared_but_not_referenced
#endif

namespace internal_451fd042b8d7665063ea81b98853365b {
static FUNC_PREFIX void
dynamic_ubb_double_precision_boundary_Dynamic_UBB_double_precision(
    uint8_t *RESTRICT const _data_indexVector, double *RESTRICT _data_pdfs,
    int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1,
    int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3,
    int64_t indexVectorSize) {

  const int64_t f_in_inv_dir_idx[] = {0, 2,  1,  4,  3,  6,  5,  10, 9, 8,
                                      7, 16, 15, 18, 17, 12, 11, 14, 13};

  const double weights[] = {
      0.333333333333333,  0.0555555555555556, 0.0555555555555556,
      0.0555555555555556, 0.0555555555555556, 0.0555555555555556,
      0.0555555555555556, 0.0277777777777778, 0.0277777777777778,
      0.0277777777777778, 0.0277777777777778, 0.0277777777777778,
      0.0277777777777778, 0.0277777777777778, 0.0277777777777778,
      0.0277777777777778, 0.0277777777777778, 0.0277777777777778,
      0.0277777777777778};

  const int64_t neighbour_offset_x[] = {0, 0, 0, -1, 1, 0, 0, -1, 1, -1,
                                        1, 0, 0, -1, 1, 0, 0, -1, 1};
  const int64_t neighbour_offset_y[] = {0,  1, -1, 0, 0, 0, 0,  1, 1, -1,
                                        -1, 1, -1, 0, 0, 1, -1, 0, 0};
  const int64_t neighbour_offset_z[] = {0, 0, 0, 0, 0, 1,  -1, 0,  0, 0,
                                        0, 1, 1, 1, 1, -1, -1, -1, -1};

  for (int64_t ctr_0 = 0; ctr_0 < indexVectorSize; ctr_0 += 1) {
    const int64_t x = *((int32_t *)(&_data_indexVector[40 * ctr_0]));
    const int64_t y = *((int32_t *)(&_data_indexVector[40 * ctr_0 + 4]));
    const int64_t z = *((int32_t *)(&_data_indexVector[40 * ctr_0 + 8]));
    const int64_t dir = *((int32_t *)(&_data_indexVector[40 * ctr_0 + 12]));
    const double vel0Term =
        _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_0 + _stride_pdfs_1 * y +
                   _stride_pdfs_1 + _stride_pdfs_2 * z + 8 * _stride_pdfs_3] +
        _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_0 + _stride_pdfs_1 * y +
                   _stride_pdfs_2 * z + 4 * _stride_pdfs_3] +
        _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_0 + _stride_pdfs_1 * y +
                   _stride_pdfs_2 * z + _stride_pdfs_2 + 14 * _stride_pdfs_3] +
        _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_0 + _stride_pdfs_1 * y +
                   _stride_pdfs_2 * z - _stride_pdfs_2 + 18 * _stride_pdfs_3] +
        _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_0 + _stride_pdfs_1 * y -
                   _stride_pdfs_1 + _stride_pdfs_2 * z + 10 * _stride_pdfs_3];
    const double vel1Term =
        _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y + _stride_pdfs_1 +
                   _stride_pdfs_2 * z + _stride_pdfs_2 + 11 * _stride_pdfs_3] +
        _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y + _stride_pdfs_1 +
                   _stride_pdfs_2 * z + _stride_pdfs_3] +
        _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y + _stride_pdfs_1 +
                   _stride_pdfs_2 * z - _stride_pdfs_2 + 15 * _stride_pdfs_3] +
        _data_pdfs[_stride_pdfs_0 * x - _stride_pdfs_0 + _stride_pdfs_1 * y +
                   _stride_pdfs_1 + _stride_pdfs_2 * z + 7 * _stride_pdfs_3];
    const double vel2Term =
        _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y +
                   _stride_pdfs_2 * z + _stride_pdfs_2 + 5 * _stride_pdfs_3] +
        _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y - _stride_pdfs_1 +
                   _stride_pdfs_2 * z + _stride_pdfs_2 + 12 * _stride_pdfs_3] +
        _data_pdfs[_stride_pdfs_0 * x - _stride_pdfs_0 + _stride_pdfs_1 * y +
                   _stride_pdfs_2 * z + _stride_pdfs_2 + 13 * _stride_pdfs_3];
    const double rho =
        vel0Term + vel1Term + vel2Term +
        _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y +
                   _stride_pdfs_2 * z - _stride_pdfs_2 + 6 * _stride_pdfs_3] +
        _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y +
                   _stride_pdfs_2 * z] +
        _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y - _stride_pdfs_1 +
                   _stride_pdfs_2 * z + 2 * _stride_pdfs_3] +
        _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y - _stride_pdfs_1 +
                   _stride_pdfs_2 * z - _stride_pdfs_2 + 16 * _stride_pdfs_3] +
        _data_pdfs[_stride_pdfs_0 * x - _stride_pdfs_0 + _stride_pdfs_1 * y +
                   _stride_pdfs_2 * z + 3 * _stride_pdfs_3] +
        _data_pdfs[_stride_pdfs_0 * x - _stride_pdfs_0 + _stride_pdfs_1 * y +
                   _stride_pdfs_2 * z - _stride_pdfs_2 + 17 * _stride_pdfs_3] +
        _data_pdfs[_stride_pdfs_0 * x - _stride_pdfs_0 + _stride_pdfs_1 * y -
                   _stride_pdfs_1 + _stride_pdfs_2 * z + 9 * _stride_pdfs_3];
    _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_0 * neighbour_offset_x[dir] +
               _stride_pdfs_1 * y + _stride_pdfs_1 * neighbour_offset_y[dir] +
               _stride_pdfs_2 * z + _stride_pdfs_2 * neighbour_offset_z[dir] +
               _stride_pdfs_3 * f_in_inv_dir_idx[dir]] =
        -rho *
            (6.0 * *((double *)(&_data_indexVector[40 * ctr_0 + 16])) *
                 neighbour_offset_x[dir] +
             6.0 * *((double *)(&_data_indexVector[40 * ctr_0 + 24])) *
                 neighbour_offset_y[dir] +
             6.0 * *((double *)(&_data_indexVector[40 * ctr_0 + 32])) *
                 neighbour_offset_z[dir]) *
            weights[dir] +
        _data_pdfs[_stride_pdfs_0 * x + _stride_pdfs_1 * y +
                   _stride_pdfs_2 * z + _stride_pdfs_3 * dir];
  }
}
} // namespace internal_451fd042b8d7665063ea81b98853365b

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#ifdef __CUDACC__
#pragma pop
#endif

void Dynamic_UBB_double_precision::run_impl(IBlock *block,
                                            IndexVectors::Type type) {
  auto *indexVectors = block->getData<IndexVectors>(indexVectorID);
  int64_t indexVectorSize = int64_c(indexVectors->indexVector(type).size());
  if (indexVectorSize == 0)
    return;

  auto pointer = indexVectors->pointerCpu(type);

  uint8_t *_data_indexVector = reinterpret_cast<uint8_t *>(pointer);

  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);

  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()));
  double *RESTRICT _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
  const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  internal_451fd042b8d7665063ea81b98853365b::
      dynamic_ubb_double_precision_boundary_Dynamic_UBB_double_precision(
          _data_indexVector, _data_pdfs, _stride_pdfs_0, _stride_pdfs_1,
          _stride_pdfs_2, _stride_pdfs_3, indexVectorSize);
}

void Dynamic_UBB_double_precision::run(IBlock *block) {
  run_impl(block, IndexVectors::ALL);
}

void Dynamic_UBB_double_precision::inner(IBlock *block) {
  run_impl(block, IndexVectors::INNER);
}

void Dynamic_UBB_double_precision::outer(IBlock *block) {
  run_impl(block, IndexVectors::OUTER);
}

} // namespace lbm
} // namespace walberla
