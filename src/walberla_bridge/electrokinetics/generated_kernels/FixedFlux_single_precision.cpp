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
//! \\file FixedFlux_single_precision.cpp
//! \\author pystencils
//======================================================================================================================

#include <cmath>

#include "FixedFlux_single_precision.h"
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

namespace internal_11b127eed0b5044a0a655e8444f26034 {
static FUNC_PREFIX void
fixedflux_single_precision_boundary_FixedFlux_single_precision(
    float *RESTRICT _data_flux, uint8_t *RESTRICT const _data_indexVector,
    int64_t const _stride_flux_0, int64_t const _stride_flux_1,
    int64_t const _stride_flux_2, int64_t const _stride_flux_3,
    int64_t indexVectorSize) {
  for (int64_t ctr_0 = 0; ctr_0 < indexVectorSize; ctr_0 += 1) {
    const int64_t x = *((int32_t *)(&_data_indexVector[28 * ctr_0]));
    const int64_t y = *((int32_t *)(&_data_indexVector[28 * ctr_0 + 4]));
    const int64_t z = *((int32_t *)(&_data_indexVector[28 * ctr_0 + 8]));

    const int64_t cx[] = {0, 0, 0, -1, 1, 0, 0,  -1, 1,  -1, 1,  0, 0, -1,
                          1, 0, 0, -1, 1, 1, -1, 1,  -1, 1,  -1, 1, -1};
    const int64_t cy[] = {0, 1, -1, 0, 0, 0, 0, 1,  1,  -1, -1, 1,  -1, 0,
                          0, 1, -1, 0, 0, 1, 1, -1, -1, 1,  1,  -1, -1};
    const int64_t cz[] = {0, 0,  0,  0,  0,  1, -1, 0, 0, 0,  0,  1,  1, 1,
                          1, -1, -1, -1, -1, 1, 1,  1, 1, -1, -1, -1, -1};
    const int64_t invdir[] = {0,  2,  1,  4,  3,  6,  5,  10, 9,
                              8,  7,  16, 15, 18, 17, 12, 11, 14,
                              13, 26, 25, 24, 23, 22, 21, 20, 19};

    const int64_t dir = *((int32_t *)(&_data_indexVector[28 * ctr_0 + 12]));
    if (((dir) == (26))) {
      _data_flux[_stride_flux_0 * x + _stride_flux_1 * y + _stride_flux_2 * z +
                 9 * _stride_flux_3] =
          -0.11111111111111111f *
              *((float *)(&_data_indexVector[28 * ctr_0 + 16])) -
          0.11111111111111111f *
              *((float *)(&_data_indexVector[28 * ctr_0 + 20])) -
          0.11111111111111111f *
              *((float *)(&_data_indexVector[28 * ctr_0 + 24]));
    } else {
      if (((dir) == (25))) {
        _data_flux[_stride_flux_0 * x + _stride_flux_0 + _stride_flux_1 * y -
                   _stride_flux_1 + _stride_flux_2 * z - _stride_flux_2 +
                   12 * _stride_flux_3] =
            -0.11111111111111111f *
                *((float *)(&_data_indexVector[28 * ctr_0 + 16])) +
            0.11111111111111111f *
                *((float *)(&_data_indexVector[28 * ctr_0 + 20])) +
            0.11111111111111111f *
                *((float *)(&_data_indexVector[28 * ctr_0 + 24]));
      } else {
        if (((dir) == (24))) {
          _data_flux[_stride_flux_0 * x + _stride_flux_1 * y +
                     _stride_flux_2 * z + 11 * _stride_flux_3] =
              -0.11111111111111111f *
                  *((float *)(&_data_indexVector[28 * ctr_0 + 16])) -
              0.11111111111111111f *
                  *((float *)(&_data_indexVector[28 * ctr_0 + 24])) +
              0.11111111111111111f *
                  *((float *)(&_data_indexVector[28 * ctr_0 + 20]));
        } else {
          if (((dir) == (23))) {
            _data_flux[_stride_flux_0 * x + _stride_flux_0 +
                       _stride_flux_1 * y + _stride_flux_1 +
                       _stride_flux_2 * z - _stride_flux_2 +
                       10 * _stride_flux_3] =
                -0.11111111111111111f *
                    *((float *)(&_data_indexVector[28 * ctr_0 + 16])) -
                0.11111111111111111f *
                    *((float *)(&_data_indexVector[28 * ctr_0 + 20])) +
                0.11111111111111111f *
                    *((float *)(&_data_indexVector[28 * ctr_0 + 24]));
          } else {
            if (((dir) == (22))) {
              _data_flux[_stride_flux_0 * x + _stride_flux_1 * y +
                         _stride_flux_2 * z + 10 * _stride_flux_3] =
                  -0.11111111111111111f *
                      *((float *)(&_data_indexVector[28 * ctr_0 + 16])) -
                  0.11111111111111111f *
                      *((float *)(&_data_indexVector[28 * ctr_0 + 20])) +
                  0.11111111111111111f *
                      *((float *)(&_data_indexVector[28 * ctr_0 + 24]));
            } else {
              if (((dir) == (21))) {
                _data_flux[_stride_flux_0 * x + _stride_flux_0 +
                           _stride_flux_1 * y - _stride_flux_1 +
                           _stride_flux_2 * z + _stride_flux_2 +
                           11 * _stride_flux_3] =
                    -0.11111111111111111f *
                        *((float *)(&_data_indexVector[28 * ctr_0 + 16])) -
                    0.11111111111111111f *
                        *((float *)(&_data_indexVector[28 * ctr_0 + 24])) +
                    0.11111111111111111f *
                        *((float *)(&_data_indexVector[28 * ctr_0 + 20]));
              } else {
                if (((dir) == (20))) {
                  _data_flux[_stride_flux_0 * x + _stride_flux_1 * y +
                             _stride_flux_2 * z + 12 * _stride_flux_3] =
                      -0.11111111111111111f *
                          *((float *)(&_data_indexVector[28 * ctr_0 + 16])) +
                      0.11111111111111111f *
                          *((float *)(&_data_indexVector[28 * ctr_0 + 20])) +
                      0.11111111111111111f *
                          *((float *)(&_data_indexVector[28 * ctr_0 + 24]));
                } else {
                  if (((dir) == (19))) {
                    _data_flux[_stride_flux_0 * x + _stride_flux_0 +
                               _stride_flux_1 * y + _stride_flux_1 +
                               _stride_flux_2 * z + _stride_flux_2 +
                               9 * _stride_flux_3] =
                        -0.11111111111111111f *
                            *((float *)(&_data_indexVector[28 * ctr_0 + 16])) -
                        0.11111111111111111f *
                            *((float *)(&_data_indexVector[28 * ctr_0 + 20])) -
                        0.11111111111111111f *
                            *((float *)(&_data_indexVector[28 * ctr_0 + 24]));
                  } else {
                    if (((dir) == (18))) {
                      _data_flux[_stride_flux_0 * x + _stride_flux_0 +
                                 _stride_flux_1 * y + _stride_flux_2 * z -
                                 _stride_flux_2 + 6 * _stride_flux_3] =
                          -0.11111111111111111f *
                              *((float
                                     *)(&_data_indexVector[28 * ctr_0 + 16])) +
                          0.11111111111111111f *
                              *((float *)(&_data_indexVector[28 * ctr_0 + 24]));
                    } else {
                      if (((dir) == (17))) {
                        _data_flux[_stride_flux_0 * x + _stride_flux_1 * y +
                                   _stride_flux_2 * z + 5 * _stride_flux_3] =
                            -0.11111111111111111f *
                                *((float *)(&_data_indexVector[28 * ctr_0 +
                                                               16])) -
                            0.11111111111111111f *
                                *((float
                                       *)(&_data_indexVector[28 * ctr_0 + 24]));
                      } else {
                        if (((dir) == (16))) {
                          _data_flux[_stride_flux_0 * x + _stride_flux_1 * y +
                                     _stride_flux_2 * z + 7 * _stride_flux_3] =
                              -0.11111111111111111f *
                                  *((float *)(&_data_indexVector[28 * ctr_0 +
                                                                 20])) -
                              0.11111111111111111f *
                                  *((float *)(&_data_indexVector[28 * ctr_0 +
                                                                 24]));
                        } else {
                          if (((dir) == (15))) {
                            _data_flux[_stride_flux_0 * x + _stride_flux_1 * y +
                                       _stride_flux_1 + _stride_flux_2 * z -
                                       _stride_flux_2 + 8 * _stride_flux_3] =
                                -0.11111111111111111f *
                                    *((float *)(&_data_indexVector[28 * ctr_0 +
                                                                   20])) +
                                0.11111111111111111f *
                                    *((float *)(&_data_indexVector[28 * ctr_0 +
                                                                   24]));
                          } else {
                            if (((dir) == (14))) {
                              _data_flux[_stride_flux_0 * x + _stride_flux_0 +
                                         _stride_flux_1 * y +
                                         _stride_flux_2 * z + _stride_flux_2 +
                                         5 * _stride_flux_3] =
                                  -0.11111111111111111f *
                                      *((float
                                             *)(&_data_indexVector[28 * ctr_0 +
                                                                   16])) -
                                  0.11111111111111111f *
                                      *((float
                                             *)(&_data_indexVector[28 * ctr_0 +
                                                                   24]));
                            } else {
                              if (((dir) == (13))) {
                                _data_flux[_stride_flux_0 * x +
                                           _stride_flux_1 * y +
                                           _stride_flux_2 * z +
                                           6 * _stride_flux_3] =
                                    -0.11111111111111111f *
                                        *((float *)(&_data_indexVector
                                                        [28 * ctr_0 + 16])) +
                                    0.11111111111111111f *
                                        *((float *)(&_data_indexVector
                                                        [28 * ctr_0 + 24]));
                              } else {
                                if (((dir) == (12))) {
                                  _data_flux[_stride_flux_0 * x +
                                             _stride_flux_1 * y +
                                             _stride_flux_2 * z +
                                             8 * _stride_flux_3] =
                                      -0.11111111111111111f *
                                          *((float *)(&_data_indexVector
                                                          [28 * ctr_0 + 20])) +
                                      0.11111111111111111f *
                                          *((float *)(&_data_indexVector
                                                          [28 * ctr_0 + 24]));
                                } else {
                                  if (((dir) == (11))) {
                                    _data_flux[_stride_flux_0 * x +
                                               _stride_flux_1 * y +
                                               _stride_flux_1 +
                                               _stride_flux_2 * z +
                                               _stride_flux_2 +
                                               7 * _stride_flux_3] =
                                        -0.11111111111111111f *
                                            *((float
                                                   *)(&_data_indexVector
                                                          [28 * ctr_0 + 20])) -
                                        0.11111111111111111f *
                                            *((float *)(&_data_indexVector
                                                            [28 * ctr_0 + 24]));
                                  } else {
                                    if (((dir) == (10))) {
                                      _data_flux[_stride_flux_0 * x +
                                                 _stride_flux_0 +
                                                 _stride_flux_1 * y -
                                                 _stride_flux_1 +
                                                 _stride_flux_2 * z +
                                                 4 * _stride_flux_3] =
                                          -0.11111111111111111f *
                                              *((float *)(&_data_indexVector
                                                              [28 * ctr_0 +
                                                               16])) +
                                          0.11111111111111111f *
                                              *((float
                                                     *)(&_data_indexVector
                                                            [28 * ctr_0 + 20]));
                                    } else {
                                      if (((dir) == (9))) {
                                        _data_flux[_stride_flux_0 * x +
                                                   _stride_flux_1 * y +
                                                   _stride_flux_2 * z +
                                                   3 * _stride_flux_3] =
                                            -0.11111111111111111f *
                                                *((float *)(&_data_indexVector
                                                                [28 * ctr_0 +
                                                                 16])) -
                                            0.11111111111111111f *
                                                *((float *)(&_data_indexVector
                                                                [28 * ctr_0 +
                                                                 20]));
                                      } else {
                                        if (((dir) == (8))) {
                                          _data_flux[_stride_flux_0 * x +
                                                     _stride_flux_0 +
                                                     _stride_flux_1 * y +
                                                     _stride_flux_1 +
                                                     _stride_flux_2 * z +
                                                     3 * _stride_flux_3] =
                                              -0.11111111111111111f *
                                                  *((float *)(&_data_indexVector
                                                                  [28 * ctr_0 +
                                                                   16])) -
                                              0.11111111111111111f *
                                                  *((float *)(&_data_indexVector
                                                                  [28 * ctr_0 +
                                                                   20]));
                                        } else {
                                          if (((dir) == (7))) {
                                            _data_flux[_stride_flux_0 * x +
                                                       _stride_flux_1 * y +
                                                       _stride_flux_2 * z +
                                                       4 * _stride_flux_3] =
                                                -0.11111111111111111f *
                                                    *((float
                                                           *)(&_data_indexVector
                                                                  [28 * ctr_0 +
                                                                   16])) +
                                                0.11111111111111111f *
                                                    *((float
                                                           *)(&_data_indexVector
                                                                  [28 * ctr_0 +
                                                                   20]));
                                          } else {
                                            if (((dir) == (6))) {
                                              _data_flux[_stride_flux_0 * x +
                                                         _stride_flux_1 * y +
                                                         _stride_flux_2 * z +
                                                         2 * _stride_flux_3] =
                                                  -0.11111111111111111f *
                                                  *((float *)(&_data_indexVector
                                                                  [28 * ctr_0 +
                                                                   24]));
                                            } else {
                                              if (((dir) == (5))) {
                                                _data_flux[_stride_flux_0 * x +
                                                           _stride_flux_1 * y +
                                                           _stride_flux_2 * z +
                                                           _stride_flux_2 +
                                                           2 * _stride_flux_3] =
                                                    -0.11111111111111111f *
                                                    *((float
                                                           *)(&_data_indexVector
                                                                  [28 * ctr_0 +
                                                                   24]));
                                              } else {
                                                if (((dir) == (4))) {
                                                  _data_flux[_stride_flux_0 *
                                                                 x +
                                                             _stride_flux_0 +
                                                             _stride_flux_1 *
                                                                 y +
                                                             _stride_flux_2 *
                                                                 z] =
                                                      -0.11111111111111111f *
                                                      *((float
                                                             *)(&_data_indexVector
                                                                    [28 *
                                                                         ctr_0 +
                                                                     16]));
                                                } else {
                                                  if (((dir) == (3))) {
                                                    _data_flux[_stride_flux_0 *
                                                                   x +
                                                               _stride_flux_1 *
                                                                   y +
                                                               _stride_flux_2 *
                                                                   z] =
                                                        -0.11111111111111111f *
                                                        *((float
                                                               *)(&_data_indexVector
                                                                      [28 *
                                                                           ctr_0 +
                                                                       16]));
                                                  } else {
                                                    if (((dir) == (2))) {
                                                      _data_flux[_stride_flux_0 *
                                                                     x +
                                                                 _stride_flux_1 *
                                                                     y +
                                                                 _stride_flux_2 *
                                                                     z +
                                                                 _stride_flux_3] =
                                                          -0.11111111111111111f *
                                                          *((float
                                                                 *)(&_data_indexVector
                                                                        [28 *
                                                                             ctr_0 +
                                                                         20]));
                                                    } else {
                                                      if (((dir) == (1))) {
                                                        _data_flux[_stride_flux_0 *
                                                                       x +
                                                                   _stride_flux_1 *
                                                                       y +
                                                                   _stride_flux_1 +
                                                                   _stride_flux_2 *
                                                                       z +
                                                                   _stride_flux_3] =
                                                            -0.11111111111111111f *
                                                            *((float
                                                                   *)(&_data_indexVector
                                                                          [28 *
                                                                               ctr_0 +
                                                                           20]));
                                                      }
                                                    }
                                                  }
                                                }
                                              }
                                            }
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
} // namespace internal_11b127eed0b5044a0a655e8444f26034

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#ifdef __CUDACC__
#pragma pop
#endif

void FixedFlux_single_precision::run_impl(IBlock *block,
                                          IndexVectors::Type type) {
  auto *indexVectors = block->uncheckedFastGetData<IndexVectors>(indexVectorID);
  int64_t indexVectorSize = int64_c(indexVectors->indexVector(type).size());
  if (indexVectorSize == 0)
    return;

  auto pointer = indexVectors->pointerCpu(type);

  uint8_t *_data_indexVector = reinterpret_cast<uint8_t *>(pointer);

  auto flux = block->getData<field::GhostLayerField<float, 13>>(fluxID);

  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(flux->nrOfGhostLayers()));
  float *RESTRICT _data_flux = flux->dataAt(0, 0, 0, 0);
  const int64_t _stride_flux_0 = int64_t(flux->xStride());
  const int64_t _stride_flux_1 = int64_t(flux->yStride());
  const int64_t _stride_flux_2 = int64_t(flux->zStride());
  const int64_t _stride_flux_3 = int64_t(1 * int64_t(flux->fStride()));
  internal_11b127eed0b5044a0a655e8444f26034::
      fixedflux_single_precision_boundary_FixedFlux_single_precision(
          _data_flux, _data_indexVector, _stride_flux_0, _stride_flux_1,
          _stride_flux_2, _stride_flux_3, indexVectorSize);
}

void FixedFlux_single_precision::run(IBlock *block) {
  run_impl(block, IndexVectors::ALL);
}

void FixedFlux_single_precision::inner(IBlock *block) {
  run_impl(block, IndexVectors::INNER);
}

void FixedFlux_single_precision::outer(IBlock *block) {
  run_impl(block, IndexVectors::OUTER);
}

} // namespace pystencils
} // namespace walberla
