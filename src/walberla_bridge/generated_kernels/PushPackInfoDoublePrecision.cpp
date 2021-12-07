// kernel generated with pystencils v0.4.3+12.g29e0e84, lbmpy v0.4.3+2.g0e17e61,
// lbmpy_walberla/pystencils_walberla from commit
// 55e6cf598e7e55f496ffaecd40bde632de76930e

#include "PushPackInfoDoublePrecision.h"
#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "stencil/Directions.h"

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

namespace walberla {
namespace lbm {

using walberla::cell::CellInterval;
using walberla::stencil::Direction;

namespace internal_pack_SW {
static FUNC_PREFIX void
pack_SW(double *RESTRICT _data_buffer, double *RESTRICT const _data_pdfs,
        int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
        int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
        int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
        int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_20_39_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_20_39;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     ((ctr_0) / (1))] =
            _data_pdfs_20_39_1m1[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
      }
    }
  }
}
} // namespace internal_pack_SW

namespace internal_pack_BW {
static FUNC_PREFIX void
pack_BW(double *RESTRICT _data_buffer, double *RESTRICT const _data_pdfs,
        int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
        int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
        int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
        int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_2m1_317 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 17 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_2m1_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_2m1_317;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     ((ctr_0) / (1))] =
            _data_pdfs_2m1_317_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
      }
    }
  }
}
} // namespace internal_pack_BW

namespace internal_pack_W {
static FUNC_PREFIX void
pack_W(double *RESTRICT _data_buffer, double *RESTRICT const _data_pdfs,
       int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
       int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
       int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
       int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_2m1_317 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 17 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_33 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_21_313 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 13 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_20_39_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_20_39;
      double *RESTRICT _data_pdfs_2m1_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_2m1_317;
      double *RESTRICT _data_pdfs_20_33_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      double *RESTRICT _data_pdfs_21_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_21_313;
      double *RESTRICT _data_pdfs_20_37_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_20_37;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1))] =
            _data_pdfs_20_39_1m1[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1)) + 1] =
            _data_pdfs_2m1_317_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1)) + 2] =
            _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1)) + 3] =
            _data_pdfs_21_313_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1)) + 4] =
            _data_pdfs_20_37_11[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
      }
    }
  }
}
} // namespace internal_pack_W

namespace internal_pack_TW {
static FUNC_PREFIX void
pack_TW(double *RESTRICT _data_buffer, double *RESTRICT const _data_pdfs,
        int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
        int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
        int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
        int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_21_313 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 13 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_21_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_21_313;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     ((ctr_0) / (1))] =
            _data_pdfs_21_313_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
      }
    }
  }
}
} // namespace internal_pack_TW

namespace internal_pack_NW {
static FUNC_PREFIX void
pack_NW(double *RESTRICT _data_buffer, double *RESTRICT const _data_pdfs,
        int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
        int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
        int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
        int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_20_37_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_20_37;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     ((ctr_0) / (1))] =
            _data_pdfs_20_37_11[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
      }
    }
  }
}
} // namespace internal_pack_NW

namespace internal_pack_BS {
static FUNC_PREFIX void
pack_BS(double *RESTRICT _data_buffer, double *RESTRICT const _data_pdfs,
        int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
        int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
        int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
        int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_2m1_316 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 16 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_2m1_316_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_2m1_316;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     ((ctr_0) / (1))] =
            _data_pdfs_2m1_316_1m1[_stride_pdfs_0 * ctr_0];
      }
    }
  }
}
} // namespace internal_pack_BS

namespace internal_pack_S {
static FUNC_PREFIX void
pack_S(double *RESTRICT _data_buffer, double *RESTRICT const _data_pdfs,
       int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
       int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
       int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
       int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_2m1_316 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 16 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_32 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_21_312 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 12 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_20_39_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_20_39;
      double *RESTRICT _data_pdfs_2m1_316_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_2m1_316;
      double *RESTRICT _data_pdfs_20_32_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_20_32;
      double *RESTRICT _data_pdfs_21_312_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_21_312;
      double *RESTRICT _data_pdfs_20_310_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_20_310;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1))] =
            _data_pdfs_20_39_1m1[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1)) + 1] =
            _data_pdfs_2m1_316_1m1[_stride_pdfs_0 * ctr_0];
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1)) + 2] =
            _data_pdfs_20_32_1m1[_stride_pdfs_0 * ctr_0];
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1)) + 3] =
            _data_pdfs_21_312_1m1[_stride_pdfs_0 * ctr_0];
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1)) + 4] =
            _data_pdfs_20_310_1m1[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
      }
    }
  }
}
} // namespace internal_pack_S

namespace internal_pack_TS {
static FUNC_PREFIX void
pack_TS(double *RESTRICT _data_buffer, double *RESTRICT const _data_pdfs,
        int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
        int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
        int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
        int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_21_312 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 12 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_21_312_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_21_312;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     ((ctr_0) / (1))] =
            _data_pdfs_21_312_1m1[_stride_pdfs_0 * ctr_0];
      }
    }
  }
}
} // namespace internal_pack_TS

namespace internal_pack_B {
static FUNC_PREFIX void
pack_B(double *RESTRICT _data_buffer, double *RESTRICT const _data_pdfs,
       int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
       int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
       int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
       int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_2m1_317 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 17 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_2m1_316 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 16 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_2m1_36 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                         _stride_pdfs_2 + 6 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_2m1_315 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 15 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_2m1_318 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 18 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_2m1_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_2m1_317;
      double *RESTRICT _data_pdfs_2m1_316_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_2m1_316;
      double *RESTRICT _data_pdfs_2m1_36_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_2m1_36;
      double *RESTRICT _data_pdfs_2m1_315_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_2m1_315;
      double *RESTRICT _data_pdfs_2m1_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_2m1_318;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1))] =
            _data_pdfs_2m1_317_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1)) + 1] =
            _data_pdfs_2m1_316_1m1[_stride_pdfs_0 * ctr_0];
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1)) + 2] =
            _data_pdfs_2m1_36_10[_stride_pdfs_0 * ctr_0];
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1)) + 3] =
            _data_pdfs_2m1_315_11[_stride_pdfs_0 * ctr_0];
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1)) + 4] =
            _data_pdfs_2m1_318_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
      }
    }
  }
}
} // namespace internal_pack_B

namespace internal_pack_T {
static FUNC_PREFIX void
pack_T(double *RESTRICT _data_buffer, double *RESTRICT const _data_pdfs,
       int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
       int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
       int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
       int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_21_313 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 13 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_21_312 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 12 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_21_35 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                        _stride_pdfs_2 + 5 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_21_311 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 11 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_21_314 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 14 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_21_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_21_313;
      double *RESTRICT _data_pdfs_21_312_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_21_312;
      double *RESTRICT _data_pdfs_21_35_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_21_35;
      double *RESTRICT _data_pdfs_21_311_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_21_311;
      double *RESTRICT _data_pdfs_21_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_21_314;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1))] =
            _data_pdfs_21_313_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1)) + 1] =
            _data_pdfs_21_312_1m1[_stride_pdfs_0 * ctr_0];
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1)) + 2] =
            _data_pdfs_21_35_10[_stride_pdfs_0 * ctr_0];
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1)) + 3] =
            _data_pdfs_21_311_11[_stride_pdfs_0 * ctr_0];
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1)) + 4] =
            _data_pdfs_21_314_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
      }
    }
  }
}
} // namespace internal_pack_T

namespace internal_pack_BN {
static FUNC_PREFIX void
pack_BN(double *RESTRICT _data_buffer, double *RESTRICT const _data_pdfs,
        int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
        int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
        int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
        int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_2m1_315 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 15 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_2m1_315_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_2m1_315;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     ((ctr_0) / (1))] =
            _data_pdfs_2m1_315_11[_stride_pdfs_0 * ctr_0];
      }
    }
  }
}
} // namespace internal_pack_BN

namespace internal_pack_N {
static FUNC_PREFIX void
pack_N(double *RESTRICT _data_buffer, double *RESTRICT const _data_pdfs,
       int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
       int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
       int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
       int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_2m1_315 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 15 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_31 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    double *RESTRICT _data_pdfs_21_311 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 11 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_20_37_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_20_37;
      double *RESTRICT _data_pdfs_2m1_315_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_2m1_315;
      double *RESTRICT _data_pdfs_20_31_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_20_31;
      double *RESTRICT _data_pdfs_21_311_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_21_311;
      double *RESTRICT _data_pdfs_20_38_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_20_38;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1))] =
            _data_pdfs_20_37_11[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0];
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1)) + 1] =
            _data_pdfs_2m1_315_11[_stride_pdfs_0 * ctr_0];
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1)) + 2] =
            _data_pdfs_20_31_11[_stride_pdfs_0 * ctr_0];
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1)) + 3] =
            _data_pdfs_21_311_11[_stride_pdfs_0 * ctr_0];
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1)) + 4] =
            _data_pdfs_20_38_11[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
      }
    }
  }
}
} // namespace internal_pack_N

namespace internal_pack_TN {
static FUNC_PREFIX void
pack_TN(double *RESTRICT _data_buffer, double *RESTRICT const _data_pdfs,
        int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
        int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
        int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
        int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_21_311 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 11 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_21_311_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_21_311;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     ((ctr_0) / (1))] =
            _data_pdfs_21_311_11[_stride_pdfs_0 * ctr_0];
      }
    }
  }
}
} // namespace internal_pack_TN

namespace internal_pack_SE {
static FUNC_PREFIX void
pack_SE(double *RESTRICT _data_buffer, double *RESTRICT const _data_pdfs,
        int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
        int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
        int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
        int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_20_310_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_20_310;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     ((ctr_0) / (1))] =
            _data_pdfs_20_310_1m1[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
      }
    }
  }
}
} // namespace internal_pack_SE

namespace internal_pack_BE {
static FUNC_PREFIX void
pack_BE(double *RESTRICT _data_buffer, double *RESTRICT const _data_pdfs,
        int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
        int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
        int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
        int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_2m1_318 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 18 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_2m1_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_2m1_318;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     ((ctr_0) / (1))] =
            _data_pdfs_2m1_318_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
      }
    }
  }
}
} // namespace internal_pack_BE

namespace internal_pack_E {
static FUNC_PREFIX void
pack_E(double *RESTRICT _data_buffer, double *RESTRICT const _data_pdfs,
       int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
       int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
       int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
       int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_2m1_318 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 18 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_34 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_21_314 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 14 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_20_310_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_20_310;
      double *RESTRICT _data_pdfs_2m1_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_2m1_318;
      double *RESTRICT _data_pdfs_20_34_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      double *RESTRICT _data_pdfs_21_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_21_314;
      double *RESTRICT _data_pdfs_20_38_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_20_38;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1))] =
            _data_pdfs_20_310_1m1[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1)) + 1] =
            _data_pdfs_2m1_318_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1)) + 2] =
            _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1)) + 3] =
            _data_pdfs_21_314_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
        _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     5 * ((ctr_0) / (1)) + 4] =
            _data_pdfs_20_38_11[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
      }
    }
  }
}
} // namespace internal_pack_E

namespace internal_pack_TE {
static FUNC_PREFIX void
pack_TE(double *RESTRICT _data_buffer, double *RESTRICT const _data_pdfs,
        int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
        int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
        int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
        int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_21_314 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 14 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_21_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_21_314;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     ((ctr_0) / (1))] =
            _data_pdfs_21_314_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
      }
    }
  }
}
} // namespace internal_pack_TE

namespace internal_pack_NE {
static FUNC_PREFIX void
pack_NE(double *RESTRICT _data_buffer, double *RESTRICT const _data_pdfs,
        int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
        int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
        int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
        int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_20_38_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_20_38;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_buffer[((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                         ((ctr_2) / (1)) +
                     ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                     ((ctr_0) / (1))] =
            _data_pdfs_20_38_11[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0];
      }
    }
  }
}
} // namespace internal_pack_NE

namespace internal_unpack_SW {
static FUNC_PREFIX void
unpack_SW(double *RESTRICT const _data_buffer, double *RESTRICT _data_pdfs,
          int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
          int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
          int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
          int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_20_39_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_20_39;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs_20_39_1m1[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] =
            _data_buffer[((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         ((ctr_0) / (1))];
      }
    }
  }
}
} // namespace internal_unpack_SW

namespace internal_unpack_BW {
static FUNC_PREFIX void
unpack_BW(double *RESTRICT const _data_buffer, double *RESTRICT _data_pdfs,
          int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
          int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
          int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
          int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_2m1_317 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 17 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_2m1_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_2m1_317;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs_2m1_317_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] =
            _data_buffer[((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         ((ctr_0) / (1))];
      }
    }
  }
}
} // namespace internal_unpack_BW

namespace internal_unpack_W {
static FUNC_PREFIX void
unpack_W(double *RESTRICT const _data_buffer, double *RESTRICT _data_pdfs,
         int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
         int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
         int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
         int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_2m1_317 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 17 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_33 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_21_313 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 13 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_20_39_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_20_39;
      double *RESTRICT _data_pdfs_2m1_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_2m1_317;
      double *RESTRICT _data_pdfs_20_33_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      double *RESTRICT _data_pdfs_21_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_21_313;
      double *RESTRICT _data_pdfs_20_37_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_20_37;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs_20_39_1m1[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1))];
        _data_pdfs_2m1_317_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1)) + 1];
        _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1)) + 2];
        _data_pdfs_21_313_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1)) + 3];
        _data_pdfs_20_37_11[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1)) + 4];
      }
    }
  }
}
} // namespace internal_unpack_W

namespace internal_unpack_TW {
static FUNC_PREFIX void
unpack_TW(double *RESTRICT const _data_buffer, double *RESTRICT _data_pdfs,
          int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
          int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
          int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
          int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_21_313 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 13 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_21_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_21_313;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs_21_313_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] =
            _data_buffer[((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         ((ctr_0) / (1))];
      }
    }
  }
}
} // namespace internal_unpack_TW

namespace internal_unpack_NW {
static FUNC_PREFIX void
unpack_NW(double *RESTRICT const _data_buffer, double *RESTRICT _data_pdfs,
          int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
          int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
          int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
          int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_20_37_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_20_37;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs_20_37_11[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] =
            _data_buffer[((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         ((ctr_0) / (1))];
      }
    }
  }
}
} // namespace internal_unpack_NW

namespace internal_unpack_BS {
static FUNC_PREFIX void
unpack_BS(double *RESTRICT const _data_buffer, double *RESTRICT _data_pdfs,
          int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
          int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
          int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
          int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_2m1_316 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 16 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_2m1_316_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_2m1_316;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs_2m1_316_1m1[_stride_pdfs_0 * ctr_0] =
            _data_buffer[((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         ((ctr_0) / (1))];
      }
    }
  }
}
} // namespace internal_unpack_BS

namespace internal_unpack_S {
static FUNC_PREFIX void
unpack_S(double *RESTRICT const _data_buffer, double *RESTRICT _data_pdfs,
         int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
         int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
         int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
         int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_2m1_316 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 16 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_32 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_21_312 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 12 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_20_39_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_20_39;
      double *RESTRICT _data_pdfs_2m1_316_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_2m1_316;
      double *RESTRICT _data_pdfs_20_32_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_20_32;
      double *RESTRICT _data_pdfs_21_312_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_21_312;
      double *RESTRICT _data_pdfs_20_310_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_20_310;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs_20_39_1m1[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1))];
        _data_pdfs_2m1_316_1m1[_stride_pdfs_0 * ctr_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1)) + 1];
        _data_pdfs_20_32_1m1[_stride_pdfs_0 * ctr_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1)) + 2];
        _data_pdfs_21_312_1m1[_stride_pdfs_0 * ctr_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1)) + 3];
        _data_pdfs_20_310_1m1[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1)) + 4];
      }
    }
  }
}
} // namespace internal_unpack_S

namespace internal_unpack_TS {
static FUNC_PREFIX void
unpack_TS(double *RESTRICT const _data_buffer, double *RESTRICT _data_pdfs,
          int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
          int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
          int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
          int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_21_312 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 12 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_21_312_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_21_312;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs_21_312_1m1[_stride_pdfs_0 * ctr_0] =
            _data_buffer[((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         ((ctr_0) / (1))];
      }
    }
  }
}
} // namespace internal_unpack_TS

namespace internal_unpack_B {
static FUNC_PREFIX void
unpack_B(double *RESTRICT const _data_buffer, double *RESTRICT _data_pdfs,
         int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
         int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
         int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
         int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_2m1_317 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 17 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_2m1_316 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 16 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_2m1_36 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                         _stride_pdfs_2 + 6 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_2m1_315 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 15 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_2m1_318 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 18 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_2m1_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_2m1_317;
      double *RESTRICT _data_pdfs_2m1_316_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_2m1_316;
      double *RESTRICT _data_pdfs_2m1_36_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_2m1_36;
      double *RESTRICT _data_pdfs_2m1_315_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_2m1_315;
      double *RESTRICT _data_pdfs_2m1_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_2m1_318;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs_2m1_317_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1))];
        _data_pdfs_2m1_316_1m1[_stride_pdfs_0 * ctr_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1)) + 1];
        _data_pdfs_2m1_36_10[_stride_pdfs_0 * ctr_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1)) + 2];
        _data_pdfs_2m1_315_11[_stride_pdfs_0 * ctr_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1)) + 3];
        _data_pdfs_2m1_318_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1)) + 4];
      }
    }
  }
}
} // namespace internal_unpack_B

namespace internal_unpack_T {
static FUNC_PREFIX void
unpack_T(double *RESTRICT const _data_buffer, double *RESTRICT _data_pdfs,
         int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
         int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
         int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
         int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_21_313 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 13 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_21_312 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 12 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_21_35 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                        _stride_pdfs_2 + 5 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_21_311 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 11 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_21_314 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 14 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_21_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_21_313;
      double *RESTRICT _data_pdfs_21_312_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_21_312;
      double *RESTRICT _data_pdfs_21_35_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_21_35;
      double *RESTRICT _data_pdfs_21_311_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_21_311;
      double *RESTRICT _data_pdfs_21_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_21_314;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs_21_313_10[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1))];
        _data_pdfs_21_312_1m1[_stride_pdfs_0 * ctr_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1)) + 1];
        _data_pdfs_21_35_10[_stride_pdfs_0 * ctr_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1)) + 2];
        _data_pdfs_21_311_11[_stride_pdfs_0 * ctr_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1)) + 3];
        _data_pdfs_21_314_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1)) + 4];
      }
    }
  }
}
} // namespace internal_unpack_T

namespace internal_unpack_BN {
static FUNC_PREFIX void
unpack_BN(double *RESTRICT const _data_buffer, double *RESTRICT _data_pdfs,
          int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
          int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
          int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
          int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_2m1_315 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 15 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_2m1_315_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_2m1_315;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs_2m1_315_11[_stride_pdfs_0 * ctr_0] =
            _data_buffer[((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         ((ctr_0) / (1))];
      }
    }
  }
}
} // namespace internal_unpack_BN

namespace internal_unpack_N {
static FUNC_PREFIX void
unpack_N(double *RESTRICT const _data_buffer, double *RESTRICT _data_pdfs,
         int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
         int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
         int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
         int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_2m1_315 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 15 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_31 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    double *RESTRICT _data_pdfs_21_311 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 11 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_20_37_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_20_37;
      double *RESTRICT _data_pdfs_2m1_315_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_2m1_315;
      double *RESTRICT _data_pdfs_20_31_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_20_31;
      double *RESTRICT _data_pdfs_21_311_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_21_311;
      double *RESTRICT _data_pdfs_20_38_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_20_38;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs_20_37_11[_stride_pdfs_0 * ctr_0 - _stride_pdfs_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1))];
        _data_pdfs_2m1_315_11[_stride_pdfs_0 * ctr_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1)) + 1];
        _data_pdfs_20_31_11[_stride_pdfs_0 * ctr_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1)) + 2];
        _data_pdfs_21_311_11[_stride_pdfs_0 * ctr_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1)) + 3];
        _data_pdfs_20_38_11[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1)) + 4];
      }
    }
  }
}
} // namespace internal_unpack_N

namespace internal_unpack_TN {
static FUNC_PREFIX void
unpack_TN(double *RESTRICT const _data_buffer, double *RESTRICT _data_pdfs,
          int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
          int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
          int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
          int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_21_311 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 11 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_21_311_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_21_311;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs_21_311_11[_stride_pdfs_0 * ctr_0] =
            _data_buffer[((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         ((ctr_0) / (1))];
      }
    }
  }
}
} // namespace internal_unpack_TN

namespace internal_unpack_SE {
static FUNC_PREFIX void
unpack_SE(double *RESTRICT const _data_buffer, double *RESTRICT _data_pdfs,
          int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
          int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
          int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
          int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_20_310_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_20_310;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs_20_310_1m1[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] =
            _data_buffer[((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         ((ctr_0) / (1))];
      }
    }
  }
}
} // namespace internal_unpack_SE

namespace internal_unpack_BE {
static FUNC_PREFIX void
unpack_BE(double *RESTRICT const _data_buffer, double *RESTRICT _data_pdfs,
          int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
          int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
          int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
          int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_2m1_318 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 18 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_2m1_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_2m1_318;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs_2m1_318_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] =
            _data_buffer[((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         ((ctr_0) / (1))];
      }
    }
  }
}
} // namespace internal_unpack_BE

namespace internal_unpack_E {
static FUNC_PREFIX void
unpack_E(double *RESTRICT const _data_buffer, double *RESTRICT _data_pdfs,
         int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
         int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
         int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
         int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_2m1_318 = _data_pdfs + _stride_pdfs_2 * ctr_2 -
                                          _stride_pdfs_2 + 18 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_34 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_21_314 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 14 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_20_310_1m1 =
          _stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _data_pdfs_20_310;
      double *RESTRICT _data_pdfs_2m1_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_2m1_318;
      double *RESTRICT _data_pdfs_20_34_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      double *RESTRICT _data_pdfs_21_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_21_314;
      double *RESTRICT _data_pdfs_20_38_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_20_38;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs_20_310_1m1[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1))];
        _data_pdfs_2m1_318_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1)) + 1];
        _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1)) + 2];
        _data_pdfs_21_314_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1)) + 3];
        _data_pdfs_20_38_11[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] =
            _data_buffer[5 * ((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         5 * ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         5 * ((ctr_0) / (1)) + 4];
      }
    }
  }
}
} // namespace internal_unpack_E

namespace internal_unpack_TE {
static FUNC_PREFIX void
unpack_TE(double *RESTRICT const _data_buffer, double *RESTRICT _data_pdfs,
          int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
          int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
          int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
          int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_21_314 = _data_pdfs + _stride_pdfs_2 * ctr_2 +
                                         _stride_pdfs_2 + 14 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_21_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_21_314;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs_21_314_10[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] =
            _data_buffer[((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         ((ctr_0) / (1))];
      }
    }
  }
}
} // namespace internal_unpack_TE

namespace internal_unpack_NE {
static FUNC_PREFIX void
unpack_NE(double *RESTRICT const _data_buffer, double *RESTRICT _data_pdfs,
          int64_t const _size_pdfs_0, int64_t const _size_pdfs_1,
          int64_t const _size_pdfs_2, int64_t const _stride_pdfs_0,
          int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
          int64_t const _stride_pdfs_3) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_pdfs_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_pdfs_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_20_38_11 =
          _stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _data_pdfs_20_38;
      for (int64_t ctr_0 = 0; ctr_0 < _size_pdfs_0; ctr_0 += 1) {
        _data_pdfs_20_38_11[_stride_pdfs_0 * ctr_0 + _stride_pdfs_0] =
            _data_buffer[((_size_pdfs_0) / (1)) * ((_size_pdfs_1) / (1)) *
                             ((ctr_2) / (1)) +
                         ((_size_pdfs_0) / (1)) * ((ctr_1) / (1)) +
                         ((ctr_0) / (1))];
      }
    }
  }
}
} // namespace internal_unpack_NE

void PushPackInfoDoublePrecision::pack(Direction dir,
                                       unsigned char *byte_buffer,
                                       IBlock *block) const {
  double *buffer = reinterpret_cast<double *>(byte_buffer);

  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);

  CellInterval ci;
  pdfs->getSliceBeforeGhostLayer(dir, ci, 1, false);

  switch (dir) {
  case stencil::SW: {
    double *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT const _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_SW::pack_SW(_data_buffer, _data_pdfs, _size_pdfs_0,
                              _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0,
                              _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::BW: {
    double *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT const _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_BW::pack_BW(_data_buffer, _data_pdfs, _size_pdfs_0,
                              _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0,
                              _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::W: {
    double *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT const _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_W::pack_W(_data_buffer, _data_pdfs, _size_pdfs_0,
                            _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0,
                            _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::TW: {
    double *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT const _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_TW::pack_TW(_data_buffer, _data_pdfs, _size_pdfs_0,
                              _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0,
                              _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::NW: {
    double *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT const _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_NW::pack_NW(_data_buffer, _data_pdfs, _size_pdfs_0,
                              _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0,
                              _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::BS: {
    double *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT const _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_BS::pack_BS(_data_buffer, _data_pdfs, _size_pdfs_0,
                              _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0,
                              _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::S: {
    double *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT const _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_S::pack_S(_data_buffer, _data_pdfs, _size_pdfs_0,
                            _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0,
                            _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::TS: {
    double *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT const _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_TS::pack_TS(_data_buffer, _data_pdfs, _size_pdfs_0,
                              _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0,
                              _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::B: {
    double *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT const _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_B::pack_B(_data_buffer, _data_pdfs, _size_pdfs_0,
                            _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0,
                            _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::T: {
    double *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT const _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_T::pack_T(_data_buffer, _data_pdfs, _size_pdfs_0,
                            _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0,
                            _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::BN: {
    double *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT const _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_BN::pack_BN(_data_buffer, _data_pdfs, _size_pdfs_0,
                              _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0,
                              _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::N: {
    double *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT const _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_N::pack_N(_data_buffer, _data_pdfs, _size_pdfs_0,
                            _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0,
                            _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::TN: {
    double *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT const _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_TN::pack_TN(_data_buffer, _data_pdfs, _size_pdfs_0,
                              _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0,
                              _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::SE: {
    double *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT const _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_SE::pack_SE(_data_buffer, _data_pdfs, _size_pdfs_0,
                              _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0,
                              _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::BE: {
    double *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT const _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_BE::pack_BE(_data_buffer, _data_pdfs, _size_pdfs_0,
                              _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0,
                              _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::E: {
    double *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT const _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_E::pack_E(_data_buffer, _data_pdfs, _size_pdfs_0,
                            _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0,
                            _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::TE: {
    double *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT const _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_TE::pack_TE(_data_buffer, _data_pdfs, _size_pdfs_0,
                              _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0,
                              _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::NE: {
    double *RESTRICT _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT const _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_pack_NE::pack_NE(_data_buffer, _data_pdfs, _size_pdfs_0,
                              _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0,
                              _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  default:
    WALBERLA_ASSERT(false);
  }
}

void PushPackInfoDoublePrecision::unpack(Direction dir,
                                         unsigned char *byte_buffer,
                                         IBlock *block) const {
  double *buffer = reinterpret_cast<double *>(byte_buffer);

  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);

  CellInterval ci;
  pdfs->getGhostRegion(dir, ci, 1, false);
  auto communciationDirection = stencil::inverseDir[dir];

  switch (communciationDirection) {
  case stencil::SW: {
    double *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_SW::unpack_SW(
        _data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2,
        _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::BW: {
    double *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_BW::unpack_BW(
        _data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2,
        _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::W: {
    double *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_W::unpack_W(_data_buffer, _data_pdfs, _size_pdfs_0,
                                _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0,
                                _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::TW: {
    double *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_TW::unpack_TW(
        _data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2,
        _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::NW: {
    double *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_NW::unpack_NW(
        _data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2,
        _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::BS: {
    double *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_BS::unpack_BS(
        _data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2,
        _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::S: {
    double *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_S::unpack_S(_data_buffer, _data_pdfs, _size_pdfs_0,
                                _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0,
                                _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::TS: {
    double *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_TS::unpack_TS(
        _data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2,
        _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::B: {
    double *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_B::unpack_B(_data_buffer, _data_pdfs, _size_pdfs_0,
                                _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0,
                                _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::T: {
    double *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_T::unpack_T(_data_buffer, _data_pdfs, _size_pdfs_0,
                                _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0,
                                _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::BN: {
    double *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_BN::unpack_BN(
        _data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2,
        _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::N: {
    double *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_N::unpack_N(_data_buffer, _data_pdfs, _size_pdfs_0,
                                _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0,
                                _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::TN: {
    double *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_TN::unpack_TN(
        _data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2,
        _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::SE: {
    double *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_SE::unpack_SE(
        _data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2,
        _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::BE: {
    double *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_BE::unpack_BE(
        _data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2,
        _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::E: {
    double *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_E::unpack_E(_data_buffer, _data_pdfs, _size_pdfs_0,
                                _size_pdfs_1, _size_pdfs_2, _stride_pdfs_0,
                                _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::TE: {
    double *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_TE::unpack_TE(
        _data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2,
        _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  case stencil::NE: {
    double *RESTRICT const _data_buffer = buffer;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
    double *RESTRICT _data_pdfs =
        pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->xSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.xSize()) + 0));
    const int64_t _size_pdfs_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->ySizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.ySize()) + 0));
    const int64_t _size_pdfs_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
    WALBERLA_ASSERT_GREATER_EQUAL(pdfs->zSizeWithGhostLayer(),
                                  int64_t(cell_idx_c(ci.zSize()) + 0));
    const int64_t _size_pdfs_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
    const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
    const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
    const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
    const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
    internal_unpack_NE::unpack_NE(
        _data_buffer, _data_pdfs, _size_pdfs_0, _size_pdfs_1, _size_pdfs_2,
        _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3);
    break;
  }

  default:
    WALBERLA_ASSERT(false);
  }
}

uint_t PushPackInfoDoublePrecision::size(stencil::Direction dir,
                                         const IBlock *block) const {
  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);

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
  return ci.numCells() * elementsPerCell * sizeof(double);
}

} // namespace lbm
} // namespace walberla