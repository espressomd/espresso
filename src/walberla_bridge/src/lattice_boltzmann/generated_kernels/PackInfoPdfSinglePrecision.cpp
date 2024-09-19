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
//! \\file PackInfoPdfSinglePrecision.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.3.3, lbmpy v1.3.3, lbmpy_walberla/pystencils_walberla from waLBerla commit b0842e1a493ce19ef1bbb8d2cf382fc343970a7f

#include "PackInfoPdfSinglePrecision.h"
#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "stencil/Directions.h"

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

namespace walberla {
namespace pystencils {

using walberla::cell::CellInterval;
using walberla::stencil::Direction;

namespace internal_pack_SW {
static FUNC_PREFIX void pack_SW(float *RESTRICT _data_buffer, float *RESTRICT const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[_size_pdfs_0 * _size_pdfs_1 * ctr_2 + _size_pdfs_0 * ctr_1 + ctr_0] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3];
      }
    }
  }
}
} // namespace internal_pack_SW

namespace internal_pack_BW {
static FUNC_PREFIX void pack_BW(float *RESTRICT _data_buffer, float *RESTRICT const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[_size_pdfs_0 * _size_pdfs_1 * ctr_2 + _size_pdfs_0 * ctr_1 + ctr_0] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3];
      }
    }
  }
}
} // namespace internal_pack_BW

namespace internal_pack_W {
static FUNC_PREFIX void pack_W(float *RESTRICT _data_buffer, float *RESTRICT const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3];
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 1] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3];
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 2] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3];
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 3] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3];
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 4] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3];
      }
    }
  }
}
} // namespace internal_pack_W

namespace internal_pack_TW {
static FUNC_PREFIX void pack_TW(float *RESTRICT _data_buffer, float *RESTRICT const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[_size_pdfs_0 * _size_pdfs_1 * ctr_2 + _size_pdfs_0 * ctr_1 + ctr_0] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3];
      }
    }
  }
}
} // namespace internal_pack_TW

namespace internal_pack_NW {
static FUNC_PREFIX void pack_NW(float *RESTRICT _data_buffer, float *RESTRICT const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[_size_pdfs_0 * _size_pdfs_1 * ctr_2 + _size_pdfs_0 * ctr_1 + ctr_0] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3];
      }
    }
  }
}
} // namespace internal_pack_NW

namespace internal_pack_BS {
static FUNC_PREFIX void pack_BS(float *RESTRICT _data_buffer, float *RESTRICT const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[_size_pdfs_0 * _size_pdfs_1 * ctr_2 + _size_pdfs_0 * ctr_1 + ctr_0] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3];
      }
    }
  }
}
} // namespace internal_pack_BS

namespace internal_pack_S {
static FUNC_PREFIX void pack_S(float *RESTRICT _data_buffer, float *RESTRICT const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3];
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 1] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3];
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 2] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3];
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 3] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3];
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 4] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3];
      }
    }
  }
}
} // namespace internal_pack_S

namespace internal_pack_TS {
static FUNC_PREFIX void pack_TS(float *RESTRICT _data_buffer, float *RESTRICT const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[_size_pdfs_0 * _size_pdfs_1 * ctr_2 + _size_pdfs_0 * ctr_1 + ctr_0] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3];
      }
    }
  }
}
} // namespace internal_pack_TS

namespace internal_pack_B {
static FUNC_PREFIX void pack_B(float *RESTRICT _data_buffer, float *RESTRICT const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3];
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 1] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3];
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 2] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3];
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 3] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3];
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 4] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3];
      }
    }
  }
}
} // namespace internal_pack_B

namespace internal_pack_T {
static FUNC_PREFIX void pack_T(float *RESTRICT _data_buffer, float *RESTRICT const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3];
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 1] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3];
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 2] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3];
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 3] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3];
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 4] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3];
      }
    }
  }
}
} // namespace internal_pack_T

namespace internal_pack_BN {
static FUNC_PREFIX void pack_BN(float *RESTRICT _data_buffer, float *RESTRICT const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[_size_pdfs_0 * _size_pdfs_1 * ctr_2 + _size_pdfs_0 * ctr_1 + ctr_0] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3];
      }
    }
  }
}
} // namespace internal_pack_BN

namespace internal_pack_N {
static FUNC_PREFIX void pack_N(float *RESTRICT _data_buffer, float *RESTRICT const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3];
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 1] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3];
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 2] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3];
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 3] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3];
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 4] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3];
      }
    }
  }
}
} // namespace internal_pack_N

namespace internal_pack_TN {
static FUNC_PREFIX void pack_TN(float *RESTRICT _data_buffer, float *RESTRICT const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[_size_pdfs_0 * _size_pdfs_1 * ctr_2 + _size_pdfs_0 * ctr_1 + ctr_0] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3];
      }
    }
  }
}
} // namespace internal_pack_TN

namespace internal_pack_SE {
static FUNC_PREFIX void pack_SE(float *RESTRICT _data_buffer, float *RESTRICT const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[_size_pdfs_0 * _size_pdfs_1 * ctr_2 + _size_pdfs_0 * ctr_1 + ctr_0] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3];
      }
    }
  }
}
} // namespace internal_pack_SE

namespace internal_pack_BE {
static FUNC_PREFIX void pack_BE(float *RESTRICT _data_buffer, float *RESTRICT const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[_size_pdfs_0 * _size_pdfs_1 * ctr_2 + _size_pdfs_0 * ctr_1 + ctr_0] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3];
      }
    }
  }
}
} // namespace internal_pack_BE

namespace internal_pack_E {
static FUNC_PREFIX void pack_E(float *RESTRICT _data_buffer, float *RESTRICT const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3];
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 1] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3];
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 2] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3];
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 3] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3];
        _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 4] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3];
      }
    }
  }
}
} // namespace internal_pack_E

namespace internal_pack_TE {
static FUNC_PREFIX void pack_TE(float *RESTRICT _data_buffer, float *RESTRICT const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[_size_pdfs_0 * _size_pdfs_1 * ctr_2 + _size_pdfs_0 * ctr_1 + ctr_0] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3];
      }
    }
  }
}
} // namespace internal_pack_TE

namespace internal_pack_NE {
static FUNC_PREFIX void pack_NE(float *RESTRICT _data_buffer, float *RESTRICT const _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[_size_pdfs_0 * _size_pdfs_1 * ctr_2 + _size_pdfs_0 * ctr_1 + ctr_0] = _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3];
      }
    }
  }
}
} // namespace internal_pack_NE

namespace internal_unpack_SW {
static FUNC_PREFIX void unpack_SW(float *RESTRICT const _data_buffer, float *RESTRICT _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3] = _data_buffer[_size_pdfs_0 * _size_pdfs_1 * ctr_2 + _size_pdfs_0 * ctr_1 + ctr_0];
      }
    }
  }
}
} // namespace internal_unpack_SW

namespace internal_unpack_BW {
static FUNC_PREFIX void unpack_BW(float *RESTRICT const _data_buffer, float *RESTRICT _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3] = _data_buffer[_size_pdfs_0 * _size_pdfs_1 * ctr_2 + _size_pdfs_0 * ctr_1 + ctr_0];
      }
    }
  }
}
} // namespace internal_unpack_BW

namespace internal_unpack_W {
static FUNC_PREFIX void unpack_W(float *RESTRICT const _data_buffer, float *RESTRICT _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0];
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 1];
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 2];
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 3];
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 4];
      }
    }
  }
}
} // namespace internal_unpack_W

namespace internal_unpack_TW {
static FUNC_PREFIX void unpack_TW(float *RESTRICT const _data_buffer, float *RESTRICT _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3] = _data_buffer[_size_pdfs_0 * _size_pdfs_1 * ctr_2 + _size_pdfs_0 * ctr_1 + ctr_0];
      }
    }
  }
}
} // namespace internal_unpack_TW

namespace internal_unpack_NW {
static FUNC_PREFIX void unpack_NW(float *RESTRICT const _data_buffer, float *RESTRICT _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3] = _data_buffer[_size_pdfs_0 * _size_pdfs_1 * ctr_2 + _size_pdfs_0 * ctr_1 + ctr_0];
      }
    }
  }
}
} // namespace internal_unpack_NW

namespace internal_unpack_BS {
static FUNC_PREFIX void unpack_BS(float *RESTRICT const _data_buffer, float *RESTRICT _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3] = _data_buffer[_size_pdfs_0 * _size_pdfs_1 * ctr_2 + _size_pdfs_0 * ctr_1 + ctr_0];
      }
    }
  }
}
} // namespace internal_unpack_BS

namespace internal_unpack_S {
static FUNC_PREFIX void unpack_S(float *RESTRICT const _data_buffer, float *RESTRICT _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0];
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 1];
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 2];
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 3];
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 4];
      }
    }
  }
}
} // namespace internal_unpack_S

namespace internal_unpack_TS {
static FUNC_PREFIX void unpack_TS(float *RESTRICT const _data_buffer, float *RESTRICT _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3] = _data_buffer[_size_pdfs_0 * _size_pdfs_1 * ctr_2 + _size_pdfs_0 * ctr_1 + ctr_0];
      }
    }
  }
}
} // namespace internal_unpack_TS

namespace internal_unpack_B {
static FUNC_PREFIX void unpack_B(float *RESTRICT const _data_buffer, float *RESTRICT _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0];
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 1];
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 2];
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 3];
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 4];
      }
    }
  }
}
} // namespace internal_unpack_B

namespace internal_unpack_T {
static FUNC_PREFIX void unpack_T(float *RESTRICT const _data_buffer, float *RESTRICT _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0];
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 1];
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 2];
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 3];
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 4];
      }
    }
  }
}
} // namespace internal_unpack_T

namespace internal_unpack_BN {
static FUNC_PREFIX void unpack_BN(float *RESTRICT const _data_buffer, float *RESTRICT _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3] = _data_buffer[_size_pdfs_0 * _size_pdfs_1 * ctr_2 + _size_pdfs_0 * ctr_1 + ctr_0];
      }
    }
  }
}
} // namespace internal_unpack_BN

namespace internal_unpack_N {
static FUNC_PREFIX void unpack_N(float *RESTRICT const _data_buffer, float *RESTRICT _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0];
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 1];
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 2];
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 3];
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 4];
      }
    }
  }
}
} // namespace internal_unpack_N

namespace internal_unpack_TN {
static FUNC_PREFIX void unpack_TN(float *RESTRICT const _data_buffer, float *RESTRICT _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3] = _data_buffer[_size_pdfs_0 * _size_pdfs_1 * ctr_2 + _size_pdfs_0 * ctr_1 + ctr_0];
      }
    }
  }
}
} // namespace internal_unpack_TN

namespace internal_unpack_SE {
static FUNC_PREFIX void unpack_SE(float *RESTRICT const _data_buffer, float *RESTRICT _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3] = _data_buffer[_size_pdfs_0 * _size_pdfs_1 * ctr_2 + _size_pdfs_0 * ctr_1 + ctr_0];
      }
    }
  }
}
} // namespace internal_unpack_SE

namespace internal_unpack_BE {
static FUNC_PREFIX void unpack_BE(float *RESTRICT const _data_buffer, float *RESTRICT _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3] = _data_buffer[_size_pdfs_0 * _size_pdfs_1 * ctr_2 + _size_pdfs_0 * ctr_1 + ctr_0];
      }
    }
  }
}
} // namespace internal_unpack_BE

namespace internal_unpack_E {
static FUNC_PREFIX void unpack_E(float *RESTRICT const _data_buffer, float *RESTRICT _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0];
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 1];
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 2];
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 3];
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3] = _data_buffer[5 * _size_pdfs_0 * _size_pdfs_1 * ctr_2 + 5 * _size_pdfs_0 * ctr_1 + 5 * ctr_0 + 4];
      }
    }
  }
}
} // namespace internal_unpack_E

namespace internal_unpack_TE {
static FUNC_PREFIX void unpack_TE(float *RESTRICT const _data_buffer, float *RESTRICT _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3] = _data_buffer[_size_pdfs_0 * _size_pdfs_1 * ctr_2 + _size_pdfs_0 * ctr_1 + ctr_0];
      }
    }
  }
}
} // namespace internal_unpack_TE

namespace internal_unpack_NE {
static FUNC_PREFIX void unpack_NE(float *RESTRICT const _data_buffer, float *RESTRICT _data_pdfs, int64_t const _size_pdfs_0, int64_t const _size_pdfs_1, int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs[_stride_pdfs_0 * ctr_0 + _stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3] = _data_buffer[_size_pdfs_0 * _size_pdfs_1 * ctr_2 + _size_pdfs_0 * ctr_1 + ctr_0];
      }
    }
  }
}
} // namespace internal_unpack_NE

void PackInfoPdfSinglePrecision::pack(Direction dir, unsigned char *byte_buffer, IBlock *block) const {
  float *buffer = reinterpret_cast<float *>(byte_buffer);

  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);

  CellInterval ci;
  pdfs->getSliceBeforeGhostLayer(dir, ci, 1, false);

  switch (dir) {
  case stencil::SW: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_SW::pack_SW(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::BW: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_BW::pack_BW(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::W: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_W::pack_W(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::TW: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_TW::pack_TW(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::NW: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_NW::pack_NW(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::BS: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_BS::pack_BS(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::S: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_S::pack_S(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::TS: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_TS::pack_TS(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::B: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_B::pack_B(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::T: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_T::pack_T(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::BN: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_BN::pack_BN(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::N: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_N::pack_N(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::TN: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_TN::pack_TN(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::SE: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_SE::pack_SE(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::BE: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_BE::pack_BE(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::E: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_E::pack_E(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::TE: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_TE::pack_TE(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::NE: {
    float *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_NE::pack_NE(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  default:
    WALBERLA_ASSERT(false);
  }
}

void PackInfoPdfSinglePrecision::unpack(Direction dir, unsigned char *byte_buffer, IBlock *block) const {
  float *buffer = reinterpret_cast<float *>(byte_buffer);

  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);

  CellInterval ci;
  pdfs->getGhostRegion(dir, ci, 1, false);
  auto communciationDirection = stencil::inverseDir[dir];

  switch (communciationDirection) {
  case stencil::SW: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_SW::unpack_SW(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::BW: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_BW::unpack_BW(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::W: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_W::unpack_W(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::TW: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_TW::unpack_TW(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::NW: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_NW::unpack_NW(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::BS: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_BS::unpack_BS(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::S: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_S::unpack_S(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::TS: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_TS::unpack_TS(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::B: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_B::unpack_B(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::T: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_T::unpack_T(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::BN: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_BN::unpack_BN(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::N: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_N::unpack_N(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::TN: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_TN::unpack_TN(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::SE: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_SE::unpack_SE(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::BE: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_BE::unpack_BE(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::E: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_E::unpack_E(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::TE: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_TE::unpack_TE(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::NE: {
    float *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()))
    float *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 0))
    const int64_t _size_pdfs_0 = int64_t(int64_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 0))
    const int64_t _size_pdfs_1 = int64_t(int64_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 0))
    const int64_t _size_pdfs_2 = int64_t(int64_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_NE::unpack_NE(_data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  default:
    WALBERLA_ASSERT(false);
  }
}

uint_t PackInfoPdfSinglePrecision::size(stencil::Direction dir, const IBlock *block) const {
  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);

  CellInterval ci;
  pdfs->getGhostRegion(dir, ci, 1, false);

  uint_t elementsPerCell = 0;

  switch (dir) {
  case stencil::SW:
    elementsPerCell = 1;
    break;

  case stencil::BW:
    elementsPerCell = 1;
    break;

  case stencil::W:
    elementsPerCell = 5;
    break;

  case stencil::TW:
    elementsPerCell = 1;
    break;

  case stencil::NW:
    elementsPerCell = 1;
    break;

  case stencil::BS:
    elementsPerCell = 1;
    break;

  case stencil::S:
    elementsPerCell = 5;
    break;

  case stencil::TS:
    elementsPerCell = 1;
    break;

  case stencil::B:
    elementsPerCell = 5;
    break;

  case stencil::T:
    elementsPerCell = 5;
    break;

  case stencil::BN:
    elementsPerCell = 1;
    break;

  case stencil::N:
    elementsPerCell = 5;
    break;

  case stencil::TN:
    elementsPerCell = 1;
    break;

  case stencil::SE:
    elementsPerCell = 1;
    break;

  case stencil::BE:
    elementsPerCell = 1;
    break;

  case stencil::E:
    elementsPerCell = 5;
    break;

  case stencil::TE:
    elementsPerCell = 1;
    break;

  case stencil::NE:
    elementsPerCell = 1;
    break;

  default:
    elementsPerCell = 0;
  }
  return ci.numCells() * elementsPerCell * sizeof(float);
}

} // namespace pystencils
} // namespace walberla
