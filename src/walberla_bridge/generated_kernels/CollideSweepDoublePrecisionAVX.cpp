// kernel generated with pystencils v0.4.4, lbmpy v0.4.4,
// lbmpy_walberla/pystencils_walberla from commit
// 88f85eb7a979f81d68e76009811aeed53ec3014e

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
//! \\file CollideSweepDoublePrecisionAVX.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

#include <cmath>

#include "CollideSweepDoublePrecisionAVX.h"
#include "core/DataTypes.h"
#include "core/Macros.h"

#include <immintrin.h>

#define FUNC_PREFIX

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
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

namespace internal_87ff45d0a5307288e29e018cd62db61e {
static FUNC_PREFIX void
collidesweepdoubleprecisionavx_collidesweepdoubleprecisionavx(
    double *RESTRICT const _data_force, double *RESTRICT _data_pdfs,
    int64_t const _size_force_0, int64_t const _size_force_1,
    int64_t const _size_force_2, int64_t const _stride_force_1,
    int64_t const _stride_force_2, int64_t const _stride_force_3,
    int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
    int64_t const _stride_pdfs_3, double omega_bulk, double omega_even,
    double omega_odd, double omega_shear) {
  const double xi_29 = omega_bulk * 0.50000000000000000;
  const double xi_56 = omega_shear * 0.041666666666666667;
  const double xi_61 = omega_bulk * 0.041666666666666667;
  const double xi_72 = omega_shear * 0.12500000000000000;
  const double xi_128 = omega_odd * 0.25000000000000000;
  const double xi_134 = omega_odd * 0.083333333333333333;
  const double xi_171 = omega_shear * 0.25000000000000000;
  const double xi_194 = omega_odd * 0.041666666666666667;
  const double xi_196 = omega_odd * 0.12500000000000000;
  const double rr_0 = 0.0;
  const double xi_54 = rr_0 * 0.041666666666666667;
  const double xi_140 = rr_0 * 0.16666666666666667;
  const double xi_176 = rr_0 * 0.083333333333333333;
  for (int64_t ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_20_311 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_312 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_35 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_315 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_316 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_313 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3;
    double *RESTRICT _data_force_20_31 =
        _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    double *RESTRICT _data_pdfs_20_318 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    double *RESTRICT _data_pdfs_20_31 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_32 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    double *RESTRICT _data_force_20_32 =
        _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    double *RESTRICT _data_pdfs_20_314 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_36 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_33 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_34 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_317 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3;
    double *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    for (int64_t ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_20_311_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_311;
      double *RESTRICT _data_pdfs_20_312_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_312;
      double *RESTRICT _data_pdfs_20_35_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_35;
      double *RESTRICT _data_pdfs_20_315_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_315;
      double *RESTRICT _data_pdfs_20_316_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_316;
      double *RESTRICT _data_pdfs_20_310_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_310;
      double *RESTRICT _data_pdfs_20_37_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_37;
      double *RESTRICT _data_pdfs_20_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_313;
      double *RESTRICT _data_force_20_31_10 =
          _stride_force_1 * ctr_1 + _data_force_20_31;
      double *RESTRICT _data_pdfs_20_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_318;
      double *RESTRICT _data_pdfs_20_30_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      double *RESTRICT _data_pdfs_20_31_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_31;
      double *RESTRICT _data_pdfs_20_32_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_32;
      double *RESTRICT _data_pdfs_20_38_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_38;
      double *RESTRICT _data_force_20_32_10 =
          _stride_force_1 * ctr_1 + _data_force_20_32;
      double *RESTRICT _data_pdfs_20_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_314;
      double *RESTRICT _data_pdfs_20_39_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_39;
      double *RESTRICT _data_pdfs_20_36_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_36;
      double *RESTRICT _data_pdfs_20_33_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      double *RESTRICT _data_pdfs_20_34_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      double *RESTRICT _data_pdfs_20_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_317;
      double *RESTRICT _data_force_20_30_10 =
          _stride_force_1 * ctr_1 + _data_force_20_30;
      {
        for (int64_t ctr_0 = 0; ctr_0 < (int64_t)((_size_force_0) / (4)) * (4);
             ctr_0 += 4) {
          const __m256d xi_220 = _mm256_load_pd(&_data_pdfs_20_311_10[ctr_0]);
          const __m256d xi_221 = _mm256_load_pd(&_data_pdfs_20_312_10[ctr_0]);
          const __m256d xi_222 = _mm256_load_pd(&_data_pdfs_20_35_10[ctr_0]);
          const __m256d xi_223 = _mm256_load_pd(&_data_pdfs_20_315_10[ctr_0]);
          const __m256d xi_224 = _mm256_load_pd(&_data_pdfs_20_316_10[ctr_0]);
          const __m256d xi_225 = _mm256_load_pd(&_data_pdfs_20_310_10[ctr_0]);
          const __m256d xi_226 = _mm256_load_pd(&_data_pdfs_20_37_10[ctr_0]);
          const __m256d xi_227 = _mm256_load_pd(&_data_pdfs_20_313_10[ctr_0]);
          const __m256d xi_228 = _mm256_load_pd(&_data_force_20_31_10[ctr_0]);
          const __m256d xi_229 = _mm256_load_pd(&_data_pdfs_20_318_10[ctr_0]);
          const __m256d xi_230 = _mm256_load_pd(&_data_pdfs_20_30_10[ctr_0]);
          const __m256d xi_231 = _mm256_load_pd(&_data_pdfs_20_31_10[ctr_0]);
          const __m256d xi_232 = _mm256_load_pd(&_data_pdfs_20_32_10[ctr_0]);
          const __m256d xi_233 = _mm256_load_pd(&_data_pdfs_20_38_10[ctr_0]);
          const __m256d xi_234 = _mm256_load_pd(&_data_force_20_32_10[ctr_0]);
          const __m256d xi_235 = _mm256_load_pd(&_data_pdfs_20_314_10[ctr_0]);
          const __m256d xi_236 = _mm256_load_pd(&_data_pdfs_20_39_10[ctr_0]);
          const __m256d xi_237 = _mm256_load_pd(&_data_pdfs_20_36_10[ctr_0]);
          const __m256d xi_238 = _mm256_load_pd(&_data_pdfs_20_33_10[ctr_0]);
          const __m256d xi_239 = _mm256_load_pd(&_data_pdfs_20_34_10[ctr_0]);
          const __m256d xi_240 = _mm256_load_pd(&_data_pdfs_20_317_10[ctr_0]);
          const __m256d xi_241 = _mm256_load_pd(&_data_force_20_30_10[ctr_0]);
          const __m256d xi_0 = _mm256_add_pd(xi_229, xi_235);
          const __m256d xi_1 = _mm256_add_pd(xi_0, xi_239);
          const __m256d xi_2 =
              _mm256_add_pd(_mm256_add_pd(xi_220, xi_223), xi_231);
          const __m256d xi_3 = _mm256_add_pd(xi_221, xi_222);
          const __m256d xi_4 = _mm256_add_pd(xi_236, xi_238);
          const __m256d xi_5 = _mm256_add_pd(xi_224, xi_232);
          const __m256d xi_6 = _mm256_add_pd(xi_237, xi_240);
          const __m256d xi_9 =
              _mm256_mul_pd(xi_236, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_10 = _mm256_add_pd(
              _mm256_mul_pd(xi_226, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
              xi_9);
          const __m256d xi_11 =
              _mm256_mul_pd(xi_240, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_12 =
              _mm256_mul_pd(xi_227, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_13 =
              _mm256_mul_pd(xi_238, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_14 =
              _mm256_add_pd(_mm256_add_pd(xi_11, xi_12), xi_13);
          const __m256d xi_15 =
              _mm256_mul_pd(xi_232, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_16 =
              _mm256_mul_pd(xi_225, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_17 = _mm256_add_pd(xi_15, xi_16);
          const __m256d xi_18 =
              _mm256_mul_pd(xi_224, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_19 =
              _mm256_mul_pd(xi_221, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_20 = _mm256_add_pd(xi_18, xi_19);
          const __m256d xi_21 =
              _mm256_mul_pd(xi_229, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_22 = _mm256_add_pd(xi_11, xi_21);
          const __m256d xi_23 =
              _mm256_mul_pd(xi_223, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_24 =
              _mm256_mul_pd(xi_237, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_25 = _mm256_add_pd(
              _mm256_add_pd(_mm256_add_pd(xi_18, xi_220), xi_23), xi_24);
          const __m256d xi_30 = _mm256_mul_pd(
              xi_228, _mm256_set_pd(0.16666666666666667, 0.16666666666666667,
                                    0.16666666666666667, 0.16666666666666667));
          const __m256d xi_31 = _mm256_mul_pd(
              xi_228,
              _mm256_set_pd(0.083333333333333333, 0.083333333333333333,
                            0.083333333333333333, 0.083333333333333333));
          const __m256d xi_43 = _mm256_mul_pd(
              xi_241, _mm256_set_pd(0.16666666666666667, 0.16666666666666667,
                                    0.16666666666666667, 0.16666666666666667));
          const __m256d xi_44 = _mm256_mul_pd(
              xi_241,
              _mm256_set_pd(0.083333333333333333, 0.083333333333333333,
                            0.083333333333333333, 0.083333333333333333));
          const __m256d xi_50 = _mm256_mul_pd(
              xi_234, _mm256_set_pd(0.16666666666666667, 0.16666666666666667,
                                    0.16666666666666667, 0.16666666666666667));
          const __m256d xi_51 = _mm256_mul_pd(
              xi_234,
              _mm256_set_pd(0.083333333333333333, 0.083333333333333333,
                            0.083333333333333333, 0.083333333333333333));
          const __m256d xi_68 = _mm256_mul_pd(
              xi_228, _mm256_set_pd(0.25000000000000000, 0.25000000000000000,
                                    0.25000000000000000, 0.25000000000000000));
          const __m256d xi_73 =
              _mm256_mul_pd(xi_228, _mm256_set_pd(xi_72, xi_72, xi_72, xi_72));
          const __m256d xi_106 =
              _mm256_mul_pd(xi_230, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_107 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(xi_222, _mm256_set_pd(3.0, 3.0, 3.0, 3.0)),
                  _mm256_mul_pd(xi_237, _mm256_set_pd(3.0, 3.0, 3.0, 3.0))),
              xi_106);
          const __m256d xi_108 = _mm256_mul_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_add_pd(
                                      _mm256_mul_pd(xi_220,
                                                    _mm256_set_pd(-3.0, -3.0,
                                                                  -3.0, -3.0)),
                                      _mm256_mul_pd(xi_221,
                                                    _mm256_set_pd(-3.0, -3.0,
                                                                  -3.0, -3.0))),
                                  _mm256_mul_pd(
                                      xi_223,
                                      _mm256_set_pd(-3.0, -3.0, -3.0, -3.0))),
                              _mm256_mul_pd(xi_224, _mm256_set_pd(-3.0, -3.0,
                                                                  -3.0, -3.0))),
                          _mm256_mul_pd(xi_231,
                                        _mm256_set_pd(3.0, 3.0, 3.0, 3.0))),
                      _mm256_mul_pd(xi_232, _mm256_set_pd(3.0, 3.0, 3.0, 3.0))),
                  xi_107),
              _mm256_set_pd(omega_even, omega_even, omega_even, omega_even));
          const __m256d xi_109 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(xi_220, _mm256_set_pd(2.0, 2.0, 2.0, 2.0)),
                      _mm256_mul_pd(xi_221, _mm256_set_pd(2.0, 2.0, 2.0, 2.0))),
                  _mm256_mul_pd(xi_223, _mm256_set_pd(2.0, 2.0, 2.0, 2.0))),
              _mm256_mul_pd(xi_224, _mm256_set_pd(2.0, 2.0, 2.0, 2.0)));
          const __m256d xi_110 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(xi_238, _mm256_set_pd(5.0, 5.0, 5.0, 5.0)),
                  _mm256_mul_pd(xi_239, _mm256_set_pd(5.0, 5.0, 5.0, 5.0))),
              xi_109);
          const __m256d xi_111 = _mm256_mul_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_add_pd(
                                      _mm256_add_pd(
                                          _mm256_mul_pd(
                                              xi_227,
                                              _mm256_set_pd(-5.0, -5.0, -5.0,
                                                            -5.0)),
                                          _mm256_mul_pd(
                                              xi_229,
                                              _mm256_set_pd(-5.0, -5.0, -5.0,
                                                            -5.0))),
                                      _mm256_mul_pd(xi_231,
                                                    _mm256_set_pd(-2.0, -2.0,
                                                                  -2.0, -2.0))),
                                  _mm256_mul_pd(
                                      xi_232,
                                      _mm256_set_pd(-2.0, -2.0, -2.0, -2.0))),
                              _mm256_mul_pd(xi_235, _mm256_set_pd(-5.0, -5.0,
                                                                  -5.0, -5.0))),
                          _mm256_mul_pd(xi_240,
                                        _mm256_set_pd(-5.0, -5.0, -5.0, -5.0))),
                      xi_107),
                  xi_110),
              _mm256_set_pd(omega_even, omega_even, omega_even, omega_even));
          const __m256d xi_114 =
              _mm256_mul_pd(xi_220, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_115 = _mm256_add_pd(xi_114, xi_19);
          const __m256d xi_116 =
              _mm256_mul_pd(xi_233, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_119 =
              _mm256_mul_pd(xi_235, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_120 = _mm256_add_pd(
              _mm256_add_pd(_mm256_add_pd(xi_119, xi_12), xi_16), xi_22);
          const __m256d xi_122 =
              _mm256_mul_pd(xi_227, _mm256_set_pd(2.0, 2.0, 2.0, 2.0));
          const __m256d xi_123 =
              _mm256_mul_pd(xi_235, _mm256_set_pd(2.0, 2.0, 2.0, 2.0));
          const __m256d xi_124 = _mm256_add_pd(
              _mm256_mul_pd(xi_229, _mm256_set_pd(2.0, 2.0, 2.0, 2.0)),
              _mm256_mul_pd(xi_240, _mm256_set_pd(2.0, 2.0, 2.0, 2.0)));
          const __m256d xi_125 = _mm256_mul_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_add_pd(
                                      _mm256_add_pd(
                                          _mm256_add_pd(
                                              _mm256_add_pd(
                                                  _mm256_add_pd(
                                                      _mm256_add_pd(
                                                          _mm256_add_pd(
                                                              _mm256_mul_pd(
                                                                  xi_222,
                                                                  _mm256_set_pd(
                                                                      -4.0,
                                                                      -4.0,
                                                                      -4.0,
                                                                      -4.0)),
                                                              _mm256_mul_pd(
                                                                  xi_225,
                                                                  _mm256_set_pd(
                                                                      -7.0,
                                                                      -7.0,
                                                                      -7.0,
                                                                      -7.0))),
                                                          _mm256_mul_pd(
                                                              xi_226,
                                                              _mm256_set_pd(
                                                                  -7.0, -7.0,
                                                                  -7.0, -7.0))),
                                                      _mm256_mul_pd(
                                                          xi_231,
                                                          _mm256_set_pd(
                                                              5.0, 5.0, 5.0,
                                                              5.0))),
                                                  _mm256_mul_pd(
                                                      xi_232,
                                                      _mm256_set_pd(5.0, 5.0,
                                                                    5.0, 5.0))),
                                              _mm256_mul_pd(
                                                  xi_233,
                                                  _mm256_set_pd(-7.0, -7.0,
                                                                -7.0, -7.0))),
                                          _mm256_mul_pd(
                                              xi_236,
                                              _mm256_set_pd(-7.0, -7.0, -7.0,
                                                            -7.0))),
                                      _mm256_mul_pd(xi_237,
                                                    _mm256_set_pd(-4.0, -4.0,
                                                                  -4.0, -4.0))),
                                  xi_106),
                              xi_110),
                          xi_122),
                      xi_123),
                  xi_124),
              _mm256_set_pd(omega_even, omega_even, omega_even, omega_even));
          const __m256d xi_126 = _mm256_add_pd(xi_114, xi_221);
          const __m256d xi_127 = _mm256_add_pd(
              _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_126, xi_15), xi_224),
                            xi_23),
              xi_231);
          const __m256d xi_129 = _mm256_mul_pd(
              xi_127, _mm256_set_pd(xi_128, xi_128, xi_128, xi_128));
          const __m256d xi_130 =
              _mm256_mul_pd(xi_226, _mm256_set_pd(2.0, 2.0, 2.0, 2.0));
          const __m256d xi_131 =
              _mm256_mul_pd(xi_225, _mm256_set_pd(2.0, 2.0, 2.0, 2.0));
          const __m256d xi_132 = _mm256_add_pd(
              _mm256_mul_pd(xi_233, _mm256_set_pd(-2.0, -2.0, -2.0, -2.0)),
              _mm256_mul_pd(xi_236, _mm256_set_pd(2.0, 2.0, 2.0, 2.0)));
          const __m256d xi_133 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(
                          _mm256_add_pd(
                              _mm256_mul_pd(xi_130, _mm256_set_pd(-1.0, -1.0,
                                                                  -1.0, -1.0)),
                              xi_131),
                          xi_132),
                      xi_15),
                  xi_2),
              xi_20);
          const __m256d xi_135 = _mm256_mul_pd(
              xi_133, _mm256_set_pd(xi_134, xi_134, xi_134, xi_134));
          const __m256d xi_136 =
              _mm256_mul_pd(xi_135, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_138 = _mm256_add_pd(xi_116, xi_225);
          const __m256d xi_142 = _mm256_add_pd(xi_227, xi_240);
          const __m256d xi_146 = _mm256_mul_pd(
              xi_125,
              _mm256_set_pd(-0.019841269841269841, -0.019841269841269841,
                            -0.019841269841269841, -0.019841269841269841));
          const __m256d xi_154 = _mm256_add_pd(xi_119, xi_227);
          const __m256d xi_155 = _mm256_add_pd(
              _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_13, xi_154), xi_21),
                            xi_239),
              xi_240);
          const __m256d xi_156 = _mm256_mul_pd(
              xi_155, _mm256_set_pd(xi_128, xi_128, xi_128, xi_128));
          const __m256d xi_157 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(
                          _mm256_mul_pd(xi_131,
                                        _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                          xi_1),
                      xi_130),
                  xi_132),
              xi_14);
          const __m256d xi_158 = _mm256_mul_pd(
              xi_157, _mm256_set_pd(xi_134, xi_134, xi_134, xi_134));
          const __m256d xi_160 =
              _mm256_mul_pd(xi_158, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_161 = _mm256_add_pd(xi_223, xi_224);
          const __m256d xi_162 = _mm256_add_pd(
              _mm256_add_pd(_mm256_add_pd(xi_115, xi_161), xi_222), xi_24);
          const __m256d xi_163 = _mm256_mul_pd(
              xi_162, _mm256_set_pd(xi_128, xi_128, xi_128, xi_128));
          const __m256d xi_166 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(
                          _mm256_mul_pd(xi_122,
                                        _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                          _mm256_mul_pd(xi_123,
                                        _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                      xi_124),
                  xi_25),
              xi_3);
          const __m256d xi_167 = _mm256_mul_pd(
              xi_166, _mm256_set_pd(xi_134, xi_134, xi_134, xi_134));
          const __m256d xi_168 =
              _mm256_mul_pd(xi_167, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_170 = xi_167;
          const __m256d xi_174 = _mm256_mul_pd(
              xi_125,
              _mm256_set_pd(0.013888888888888889, 0.013888888888888889,
                            0.013888888888888889, 0.013888888888888889));
          const __m256d xi_190 = _mm256_mul_pd(
              xi_111,
              _mm256_set_pd(-0.0071428571428571429, -0.0071428571428571429,
                            -0.0071428571428571429, -0.0071428571428571429));
          const __m256d xi_192 = _mm256_mul_pd(
              xi_108,
              _mm256_set_pd(0.025000000000000000, 0.025000000000000000,
                            0.025000000000000000, 0.025000000000000000));
          const __m256d xi_195 = _mm256_mul_pd(
              xi_166, _mm256_set_pd(xi_194, xi_194, xi_194, xi_194));
          const __m256d xi_197 = _mm256_mul_pd(
              xi_162, _mm256_set_pd(xi_196, xi_196, xi_196, xi_196));
          const __m256d xi_198 = _mm256_mul_pd(
              xi_125,
              _mm256_set_pd(-0.0039682539682539683, -0.0039682539682539683,
                            -0.0039682539682539683, -0.0039682539682539683));
          const __m256d xi_202 = _mm256_mul_pd(
              xi_133, _mm256_set_pd(xi_194, xi_194, xi_194, xi_194));
          const __m256d xi_203 = _mm256_mul_pd(
              xi_127, _mm256_set_pd(xi_196, xi_196, xi_196, xi_196));
          const __m256d xi_209 = _mm256_mul_pd(
              xi_111,
              _mm256_set_pd(0.017857142857142857, 0.017857142857142857,
                            0.017857142857142857, 0.017857142857142857));
          const __m256d xi_212 = _mm256_mul_pd(
              xi_155, _mm256_set_pd(xi_196, xi_196, xi_196, xi_196));
          const __m256d xi_213 = _mm256_mul_pd(
              xi_157, _mm256_set_pd(xi_194, xi_194, xi_194, xi_194));
          const __m256d xi_32 =
              _mm256_mul_pd(xi_31, _mm256_set_pd(rr_0, rr_0, rr_0, rr_0));
          const __m256d xi_45 =
              _mm256_mul_pd(xi_44, _mm256_set_pd(rr_0, rr_0, rr_0, rr_0));
          const __m256d xi_52 =
              _mm256_mul_pd(xi_51, _mm256_set_pd(rr_0, rr_0, rr_0, rr_0));
          const __m256d xi_55 =
              _mm256_mul_pd(xi_241, _mm256_set_pd(xi_54, xi_54, xi_54, xi_54));
          const __m256d xi_60 =
              _mm256_mul_pd(xi_228, _mm256_set_pd(xi_54, xi_54, xi_54, xi_54));
          const __m256d xi_82 =
              _mm256_mul_pd(xi_234, _mm256_set_pd(xi_54, xi_54, xi_54, xi_54));
          const __m256d vel0Term =
              _mm256_add_pd(_mm256_add_pd(xi_1, xi_225), xi_233);
          const __m256d vel1Term = _mm256_add_pd(xi_2, xi_226);
          const __m256d vel2Term = _mm256_add_pd(xi_227, xi_3);
          const __m256d rho = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(
                          _mm256_add_pd(_mm256_add_pd(vel0Term, vel1Term),
                                        vel2Term),
                          xi_230),
                      xi_4),
                  xi_5),
              xi_6);
          const __m256d xi_7 =
              _mm256_div_pd(_mm256_set_pd(1.0, 1.0, 1.0, 1.0), rho);
          const __m256d xi_8 = _mm256_mul_pd(
              xi_7, _mm256_set_pd(0.50000000000000000, 0.50000000000000000,
                                  0.50000000000000000, 0.50000000000000000));
          const __m256d u_0 = _mm256_add_pd(
              _mm256_mul_pd(
                  xi_7, _mm256_add_pd(_mm256_add_pd(vel0Term, xi_10), xi_14)),
              _mm256_mul_pd(xi_241, xi_8));
          const __m256d xi_26 = _mm256_mul_pd(u_0, xi_241);
          const __m256d xi_38 = _mm256_mul_pd(
              xi_26, _mm256_set_pd(0.16666666666666667, 0.16666666666666667,
                                   0.16666666666666667, 0.16666666666666667));
          const __m256d xi_39 = _mm256_mul_pd(
              xi_26, _mm256_set_pd(0.083333333333333333, 0.083333333333333333,
                                   0.083333333333333333, 0.083333333333333333));
          const __m256d xi_40 =
              _mm256_mul_pd(xi_39, _mm256_set_pd(omega_shear, omega_shear,
                                                 omega_shear, omega_shear));
          const __m256d xi_41 = _mm256_add_pd(
              _mm256_mul_pd(xi_38, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
              xi_40);
          const __m256d xi_57 = _mm256_add_pd(
              _mm256_mul_pd(
                  _mm256_mul_pd(xi_26, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                  _mm256_set_pd(xi_56, xi_56, xi_56, xi_56)),
              xi_38);
          const __m256d xi_58 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(xi_44, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                  xi_55),
              xi_57);
          const __m256d xi_62 = _mm256_mul_pd(
              _mm256_mul_pd(xi_26, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
              _mm256_set_pd(xi_61, xi_61, xi_61, xi_61));
          const __m256d xi_69 = _mm256_mul_pd(u_0, xi_68);
          const __m256d xi_74 = _mm256_mul_pd(u_0, xi_73);
          const __m256d xi_78 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(xi_55, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                  xi_44),
              xi_57);
          const __m256d xi_85 =
              _mm256_mul_pd(xi_39, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_96 = _mm256_mul_pd(u_0, xi_234);
          const __m256d xi_97 = _mm256_mul_pd(
              xi_96, _mm256_set_pd(0.25000000000000000, 0.25000000000000000,
                                   0.25000000000000000, 0.25000000000000000));
          const __m256d xi_100 =
              _mm256_mul_pd(xi_96, _mm256_set_pd(xi_72, xi_72, xi_72, xi_72));
          const __m256d xi_112 = _mm256_mul_pd(rho, (_mm256_mul_pd(u_0, u_0)));
          const __m256d xi_151 = _mm256_mul_pd(rho, u_0);
          const __m256d xi_152 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(
                          _mm256_mul_pd(vel0Term,
                                        _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                          xi_142),
                      xi_151),
                  xi_226),
              xi_4);
          const __m256d xi_153 = _mm256_mul_pd(
              xi_152, _mm256_set_pd(xi_140, xi_140, xi_140, xi_140));
          const __m256d xi_180 = _mm256_mul_pd(
              xi_152, _mm256_set_pd(xi_176, xi_176, xi_176, xi_176));
          const __m256d u_1 = _mm256_add_pd(
              _mm256_mul_pd(
                  xi_7,
                  _mm256_add_pd(
                      _mm256_add_pd(
                          _mm256_add_pd(_mm256_add_pd(vel1Term, xi_17), xi_20),
                          xi_233),
                      xi_9)),
              _mm256_mul_pd(xi_228, xi_8));
          const __m256d xi_27 = _mm256_mul_pd(u_1, xi_228);
          const __m256d xi_33 = _mm256_mul_pd(
              xi_27, _mm256_set_pd(0.16666666666666667, 0.16666666666666667,
                                   0.16666666666666667, 0.16666666666666667));
          const __m256d xi_46 = _mm256_mul_pd(
              xi_27, _mm256_set_pd(0.083333333333333333, 0.083333333333333333,
                                   0.083333333333333333, 0.083333333333333333));
          const __m256d xi_47 =
              _mm256_mul_pd(xi_46, _mm256_set_pd(omega_shear, omega_shear,
                                                 omega_shear, omega_shear));
          const __m256d xi_48 = _mm256_add_pd(
              _mm256_mul_pd(xi_33, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
              xi_47);
          const __m256d xi_63 = _mm256_mul_pd(
              _mm256_mul_pd(xi_27, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
              _mm256_set_pd(xi_61, xi_61, xi_61, xi_61));
          const __m256d xi_70 = _mm256_mul_pd(
              u_1, _mm256_set_pd(0.25000000000000000, 0.25000000000000000,
                                 0.25000000000000000, 0.25000000000000000));
          const __m256d xi_71 = _mm256_mul_pd(xi_241, xi_70);
          const __m256d xi_75 =
              _mm256_mul_pd(u_1, _mm256_set_pd(xi_72, xi_72, xi_72, xi_72));
          const __m256d xi_76 = _mm256_mul_pd(xi_241, xi_75);
          const __m256d xi_77 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(xi_69,
                                    _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                      _mm256_mul_pd(xi_71,
                                    _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                  xi_74),
              xi_76);
          const __m256d xi_79 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(xi_74,
                                    _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                      _mm256_mul_pd(xi_76,
                                    _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                  xi_69),
              xi_71);
          const __m256d xi_87 = _mm256_mul_pd(xi_234, xi_70);
          const __m256d xi_89 = _mm256_mul_pd(xi_234, xi_75);
          const __m256d xi_94 =
              _mm256_mul_pd(xi_46, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_117 = _mm256_mul_pd(rho, (_mm256_mul_pd(u_1, u_1)));
          const __m256d xi_118 =
              _mm256_add_pd(_mm256_add_pd(xi_10, xi_116), xi_117);
          const __m256d xi_137 = _mm256_mul_pd(rho, u_1);
          const __m256d xi_139 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(
                          _mm256_add_pd(
                              _mm256_mul_pd(
                                  vel1Term,
                                  _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                              xi_137),
                          xi_138),
                      xi_221),
                  xi_236),
              xi_5);
          const __m256d xi_141 = _mm256_mul_pd(
              xi_139, _mm256_set_pd(xi_140, xi_140, xi_140, xi_140));
          const __m256d xi_172 = _mm256_mul_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(_mm256_mul_pd(u_0, xi_137), xi_138),
                      xi_226),
                  xi_9),
              _mm256_set_pd(xi_171, xi_171, xi_171, xi_171));
          const __m256d xi_177 = _mm256_mul_pd(
              xi_139, _mm256_set_pd(xi_176, xi_176, xi_176, xi_176));
          const __m256d xi_178 = xi_177;
          const __m256d xi_179 = _mm256_add_pd(xi_135, xi_178);
          const __m256d xi_188 =
              _mm256_mul_pd(xi_177, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_189 = _mm256_add_pd(xi_136, xi_188);
          const __m256d xi_204 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(xi_202, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                  xi_178),
              xi_203);
          const __m256d xi_205 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(xi_203, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                  xi_188),
              xi_202);
          const __m256d u_2 = _mm256_add_pd(
              _mm256_mul_pd(
                  xi_7,
                  _mm256_add_pd(
                      _mm256_add_pd(_mm256_add_pd(vel2Term, xi_22), xi_235),
                      xi_25)),
              _mm256_mul_pd(xi_234, xi_8));
          const __m256d xi_28 = _mm256_mul_pd(u_2, xi_234);
          const __m256d xi_34 = _mm256_mul_pd(
              xi_28, _mm256_set_pd(0.16666666666666667, 0.16666666666666667,
                                   0.16666666666666667, 0.16666666666666667));
          const __m256d xi_35 = _mm256_mul_pd(
              xi_28, _mm256_set_pd(0.083333333333333333, 0.083333333333333333,
                                   0.083333333333333333, 0.083333333333333333));
          const __m256d xi_36 =
              _mm256_mul_pd(xi_35, _mm256_set_pd(omega_shear, omega_shear,
                                                 omega_shear, omega_shear));
          const __m256d xi_37 = _mm256_add_pd(
              _mm256_mul_pd(xi_34, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
              xi_36);
          const __m256d xi_42 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(xi_27, _mm256_set_pd(0.33333333333333333,
                                                         0.33333333333333333,
                                                         0.33333333333333333,
                                                         0.33333333333333333)),
                      _mm256_mul_pd(
                          _mm256_mul_pd(xi_33,
                                        _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                          _mm256_set_pd(omega_shear, omega_shear, omega_shear,
                                        omega_shear))),
                  xi_37),
              xi_41);
          const __m256d xi_49 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(xi_26, _mm256_set_pd(0.33333333333333333,
                                                         0.33333333333333333,
                                                         0.33333333333333333,
                                                         0.33333333333333333)),
                      _mm256_mul_pd(
                          _mm256_mul_pd(xi_38,
                                        _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                          _mm256_set_pd(omega_shear, omega_shear, omega_shear,
                                        omega_shear))),
                  xi_37),
              xi_48);
          const __m256d xi_53 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(xi_28, _mm256_set_pd(0.33333333333333333,
                                                         0.33333333333333333,
                                                         0.33333333333333333,
                                                         0.33333333333333333)),
                      _mm256_mul_pd(
                          _mm256_mul_pd(xi_34,
                                        _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                          _mm256_set_pd(omega_shear, omega_shear, omega_shear,
                                        omega_shear))),
                  xi_41),
              xi_48);
          const __m256d xi_59 =
              _mm256_mul_pd(xi_35, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_64 = _mm256_mul_pd(
              _mm256_mul_pd(xi_28, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
              _mm256_set_pd(xi_61, xi_61, xi_61, xi_61));
          const __m256d xi_65 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(
                          _mm256_mul_pd(
                              _mm256_mul_pd(
                                  xi_27, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                              _mm256_set_pd(xi_56, xi_56, xi_56, xi_56)),
                          xi_33),
                      xi_62),
                  xi_63),
              xi_64);
          const __m256d xi_66 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(xi_60, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                  xi_31),
              xi_65);
          const __m256d xi_67 =
              _mm256_add_pd(_mm256_add_pd(xi_36, xi_59), xi_66);
          const __m256d xi_80 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(xi_31, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                  xi_60),
              xi_65);
          const __m256d xi_81 =
              _mm256_add_pd(_mm256_add_pd(xi_36, xi_59), xi_80);
          const __m256d xi_83 = _mm256_add_pd(
              _mm256_mul_pd(
                  _mm256_mul_pd(xi_28, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                  _mm256_set_pd(xi_56, xi_56, xi_56, xi_56)),
              xi_34);
          const __m256d xi_84 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(xi_82, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                  xi_51),
              xi_83);
          const __m256d xi_86 =
              _mm256_add_pd(_mm256_add_pd(xi_40, xi_66), xi_85);
          const __m256d xi_88 = _mm256_mul_pd(u_2, xi_68);
          const __m256d xi_90 = _mm256_mul_pd(u_2, xi_73);
          const __m256d xi_91 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(xi_89,
                                    _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                      _mm256_mul_pd(xi_90,
                                    _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                  xi_87),
              xi_88);
          const __m256d xi_92 =
              _mm256_add_pd(_mm256_add_pd(xi_40, xi_80), xi_85);
          const __m256d xi_93 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(xi_87,
                                    _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                      _mm256_mul_pd(xi_88,
                                    _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                  xi_89),
              xi_90);
          const __m256d xi_95 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(_mm256_add_pd(xi_47, xi_62), xi_63), xi_64),
                  xi_84),
              xi_94);
          const __m256d xi_98 = _mm256_mul_pd(u_2, xi_241);
          const __m256d xi_99 = _mm256_mul_pd(
              xi_98, _mm256_set_pd(0.25000000000000000, 0.25000000000000000,
                                   0.25000000000000000, 0.25000000000000000));
          const __m256d xi_101 =
              _mm256_mul_pd(xi_98, _mm256_set_pd(xi_72, xi_72, xi_72, xi_72));
          const __m256d xi_102 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(xi_97,
                                    _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                      _mm256_mul_pd(xi_99,
                                    _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                  xi_100),
              xi_101);
          const __m256d xi_103 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(xi_100,
                                    _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                      _mm256_mul_pd(xi_101,
                                    _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                  xi_97),
              xi_99);
          const __m256d xi_104 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(xi_51, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                  xi_82),
              xi_83);
          const __m256d xi_105 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(_mm256_add_pd(xi_104, xi_47), xi_62),
                      xi_63),
                  xi_64),
              xi_94);
          const __m256d xi_113 = _mm256_mul_pd(rho, (_mm256_mul_pd(u_2, u_2)));
          const __m256d xi_121 = _mm256_mul_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_add_pd(_mm256_add_pd(xi_112, xi_113),
                                                xi_115),
                                  xi_118),
                              xi_120),
                          xi_18),
                      xi_23),
                  xi_230),
              _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk, omega_bulk));
          const __m256d xi_143 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(xi_113, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                  xi_222),
              xi_237);
          const __m256d xi_144 = _mm256_mul_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(
                                      xi_231,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                                  xi_0),
                              xi_118),
                          xi_142),
                      xi_143),
                  xi_17),
              _mm256_set_pd(omega_shear, omega_shear, omega_shear,
                            omega_shear));
          const __m256d xi_145 = _mm256_mul_pd(
              xi_144, _mm256_set_pd(0.12500000000000000, 0.12500000000000000,
                                    0.12500000000000000, 0.12500000000000000));
          const __m256d xi_147 = _mm256_mul_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_add_pd(
                                      _mm256_add_pd(
                                          _mm256_add_pd(
                                              _mm256_add_pd(
                                                  _mm256_add_pd(
                                                      _mm256_mul_pd(
                                                          xi_112,
                                                          _mm256_set_pd(
                                                              2.0, 2.0, 2.0,
                                                              2.0)),
                                                      _mm256_mul_pd(
                                                          xi_117,
                                                          _mm256_set_pd(
                                                              -1.0, -1.0, -1.0,
                                                              -1.0))),
                                                  _mm256_mul_pd(
                                                      xi_238, _mm256_set_pd(
                                                                  -2.0, -2.0,
                                                                  -2.0, -2.0))),
                                              _mm256_mul_pd(
                                                  xi_239,
                                                  _mm256_set_pd(-2.0, -2.0,
                                                                -2.0, -2.0))),
                                          xi_10),
                                      xi_109),
                                  xi_116),
                              xi_120),
                          xi_143),
                      xi_231),
                  xi_232),
              _mm256_set_pd(omega_shear, omega_shear, omega_shear,
                            omega_shear));
          const __m256d xi_148 = _mm256_mul_pd(
              xi_147,
              _mm256_set_pd(-0.041666666666666667, -0.041666666666666667,
                            -0.041666666666666667, -0.041666666666666667));
          const __m256d xi_149 = _mm256_add_pd(
              _mm256_mul_pd(xi_108, _mm256_set_pd(-0.050000000000000000,
                                                  -0.050000000000000000,
                                                  -0.050000000000000000,
                                                  -0.050000000000000000)),
              xi_148);
          const __m256d xi_150 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(xi_111,
                                    _mm256_set_pd(0.014285714285714286,
                                                  0.014285714285714286,
                                                  0.014285714285714286,
                                                  0.014285714285714286)),
                      xi_145),
                  xi_146),
              xi_149);
          const __m256d xi_159 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(xi_111, _mm256_set_pd(-0.035714285714285714,
                                                      -0.035714285714285714,
                                                      -0.035714285714285714,
                                                      -0.035714285714285714)),
                  _mm256_mul_pd(xi_147, _mm256_set_pd(0.083333333333333333,
                                                      0.083333333333333333,
                                                      0.083333333333333333,
                                                      0.083333333333333333))),
              xi_146);
          const __m256d xi_164 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(rho, u_2),
                                  _mm256_mul_pd(
                                      vel2Term,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                              xi_114),
                          xi_119),
                      xi_161),
                  xi_229),
              xi_6);
          const __m256d xi_165 = _mm256_mul_pd(
              xi_164, _mm256_set_pd(xi_140, xi_140, xi_140, xi_140));
          const __m256d xi_169 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_mul_pd(xi_111,
                                    _mm256_set_pd(-0.021428571428571429,
                                                  -0.021428571428571429,
                                                  -0.021428571428571429,
                                                  -0.021428571428571429)),
                      _mm256_mul_pd(xi_125,
                                    _mm256_set_pd(0.015873015873015873,
                                                  0.015873015873015873,
                                                  0.015873015873015873,
                                                  0.015873015873015873))),
                  _mm256_mul_pd(xi_145, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
              xi_149);
          const __m256d xi_173 = _mm256_mul_pd(
              xi_144,
              _mm256_set_pd(0.062500000000000000, 0.062500000000000000,
                            0.062500000000000000, 0.062500000000000000));
          const __m256d xi_175 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(xi_172, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                  xi_173),
              xi_174);
          const __m256d xi_181 = _mm256_mul_pd(
              xi_121,
              _mm256_set_pd(0.041666666666666667, 0.041666666666666667,
                            0.041666666666666667, 0.041666666666666667));
          const __m256d xi_182 = _mm256_add_pd(
              _mm256_mul_pd(xi_147, _mm256_set_pd(0.020833333333333333,
                                                  0.020833333333333333,
                                                  0.020833333333333333,
                                                  0.020833333333333333)),
              xi_181);
          const __m256d xi_183 = _mm256_add_pd(
              _mm256_mul_pd(xi_180, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
              xi_182);
          const __m256d xi_184 = _mm256_add_pd(xi_160, xi_183);
          const __m256d xi_185 =
              _mm256_add_pd(_mm256_add_pd(xi_172, xi_173), xi_174);
          const __m256d xi_186 = _mm256_add_pd(xi_180, xi_182);
          const __m256d xi_187 = _mm256_add_pd(xi_158, xi_186);
          const __m256d xi_191 = _mm256_mul_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(_mm256_mul_pd(u_2, xi_137), xi_126), xi_18),
                  xi_223),
              _mm256_set_pd(xi_171, xi_171, xi_171, xi_171));
          const __m256d xi_193 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(_mm256_add_pd(xi_148, xi_181), xi_190), xi_191),
              xi_192);
          const __m256d xi_199 = _mm256_mul_pd(
              xi_164, _mm256_set_pd(xi_176, xi_176, xi_176, xi_176));
          const __m256d xi_200 = _mm256_add_pd(xi_198, xi_199);
          const __m256d xi_201 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(xi_195, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                  xi_197),
              xi_200);
          const __m256d xi_206 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(
                          _mm256_mul_pd(xi_191,
                                        _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                          xi_148),
                      xi_181),
                  xi_190),
              xi_192);
          const __m256d xi_207 = _mm256_mul_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(_mm256_mul_pd(u_2, xi_151), xi_11), xi_154),
                  xi_229),
              _mm256_set_pd(xi_171, xi_171, xi_171, xi_171));
          const __m256d xi_208 =
              _mm256_mul_pd(xi_173, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_210 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(xi_207, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                  xi_208),
              xi_209);
          const __m256d xi_211 = _mm256_add_pd(xi_170, xi_200);
          const __m256d xi_214 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(xi_212, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                  xi_183),
              xi_213);
          const __m256d xi_215 =
              _mm256_add_pd(_mm256_add_pd(xi_207, xi_208), xi_209);
          const __m256d xi_216 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(xi_213, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                  xi_186),
              xi_212);
          const __m256d xi_217 = _mm256_add_pd(
              _mm256_mul_pd(xi_199, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
              xi_198);
          const __m256d xi_218 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(xi_197, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                  xi_195),
              xi_217);
          const __m256d xi_219 = _mm256_add_pd(xi_168, xi_217);
          const __m256d forceTerm_0 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(
                          _mm256_add_pd(
                              _mm256_mul_pd(
                                  xi_26, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                              _mm256_mul_pd(
                                  xi_26,
                                  _mm256_set_pd(xi_29, xi_29, xi_29, xi_29))),
                          _mm256_mul_pd(xi_27,
                                        _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                      _mm256_mul_pd(xi_27,
                                    _mm256_set_pd(xi_29, xi_29, xi_29, xi_29))),
                  _mm256_mul_pd(xi_28, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
              _mm256_mul_pd(xi_28, _mm256_set_pd(xi_29, xi_29, xi_29, xi_29)));
          const __m256d forceTerm_1 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(xi_32, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                  xi_30),
              xi_42);
          const __m256d forceTerm_2 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(xi_30, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                  xi_32),
              xi_42);
          const __m256d forceTerm_3 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(xi_43, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                  xi_45),
              xi_49);
          const __m256d forceTerm_4 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(xi_45, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                  xi_43),
              xi_49);
          const __m256d forceTerm_5 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(xi_52, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                  xi_50),
              xi_53);
          const __m256d forceTerm_6 = _mm256_add_pd(
              _mm256_add_pd(
                  _mm256_mul_pd(xi_50, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                  xi_52),
              xi_53);
          const __m256d forceTerm_7 =
              _mm256_add_pd(_mm256_add_pd(xi_58, xi_67), xi_77);
          const __m256d forceTerm_8 =
              _mm256_add_pd(_mm256_add_pd(xi_67, xi_78), xi_79);
          const __m256d forceTerm_9 =
              _mm256_add_pd(_mm256_add_pd(xi_58, xi_79), xi_81);
          const __m256d forceTerm_10 =
              _mm256_add_pd(_mm256_add_pd(xi_77, xi_78), xi_81);
          const __m256d forceTerm_11 =
              _mm256_add_pd(_mm256_add_pd(xi_84, xi_86), xi_91);
          const __m256d forceTerm_12 =
              _mm256_add_pd(_mm256_add_pd(xi_84, xi_92), xi_93);
          const __m256d forceTerm_13 =
              _mm256_add_pd(_mm256_add_pd(xi_102, xi_58), xi_95);
          const __m256d forceTerm_14 =
              _mm256_add_pd(_mm256_add_pd(xi_103, xi_78), xi_95);
          const __m256d forceTerm_15 =
              _mm256_add_pd(_mm256_add_pd(xi_104, xi_86), xi_93);
          const __m256d forceTerm_16 =
              _mm256_add_pd(_mm256_add_pd(xi_104, xi_91), xi_92);
          const __m256d forceTerm_17 =
              _mm256_add_pd(_mm256_add_pd(xi_103, xi_105), xi_58);
          const __m256d forceTerm_18 =
              _mm256_add_pd(_mm256_add_pd(xi_102, xi_105), xi_78);
          _mm256_store_pd(
              &_data_pdfs_20_30_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(
                                      xi_108,
                                      _mm256_set_pd(0.10000000000000000,
                                                    0.10000000000000000,
                                                    0.10000000000000000,
                                                    0.10000000000000000)),
                                  _mm256_mul_pd(
                                      xi_111,
                                      _mm256_set_pd(0.042857142857142857,
                                                    0.042857142857142857,
                                                    0.042857142857142857,
                                                    0.042857142857142857))),
                              _mm256_mul_pd(
                                  xi_121, _mm256_set_pd(-0.50000000000000000,
                                                        -0.50000000000000000,
                                                        -0.50000000000000000,
                                                        -0.50000000000000000))),
                          _mm256_mul_pd(xi_125,
                                        _mm256_set_pd(0.023809523809523810,
                                                      0.023809523809523810,
                                                      0.023809523809523810,
                                                      0.023809523809523810))),
                      forceTerm_0),
                  xi_230));
          _mm256_store_pd(
              &_data_pdfs_20_31_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(
                                      xi_129,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                                  forceTerm_1),
                              xi_136),
                          xi_141),
                      xi_150),
                  xi_231));
          _mm256_store_pd(
              &_data_pdfs_20_32_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(
                                      xi_141,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                                  forceTerm_2),
                              xi_129),
                          xi_135),
                      xi_150),
                  xi_232));
          _mm256_store_pd(
              &_data_pdfs_20_33_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(
                                      xi_153,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                                  forceTerm_3),
                              xi_156),
                          xi_158),
                      xi_159),
                  xi_238));
          _mm256_store_pd(
              &_data_pdfs_20_34_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(
                                      xi_156,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                                  forceTerm_4),
                              xi_153),
                          xi_159),
                      xi_160),
                  xi_239));
          _mm256_store_pd(
              &_data_pdfs_20_35_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(
                                      xi_163,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                                  forceTerm_5),
                              xi_165),
                          xi_168),
                      xi_169),
                  xi_222));
          _mm256_store_pd(
              &_data_pdfs_20_36_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(
                          _mm256_add_pd(
                              _mm256_add_pd(
                                  _mm256_mul_pd(
                                      xi_165,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                                  forceTerm_6),
                              xi_163),
                          xi_169),
                      xi_170),
                  xi_237));
          _mm256_store_pd(
              &_data_pdfs_20_37_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(_mm256_add_pd(forceTerm_7, xi_175), xi_179),
                      xi_184),
                  xi_226));
          _mm256_store_pd(
              &_data_pdfs_20_38_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(_mm256_add_pd(forceTerm_8, xi_179), xi_185),
                      xi_187),
                  xi_233));
          _mm256_store_pd(
              &_data_pdfs_20_39_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(_mm256_add_pd(forceTerm_9, xi_184), xi_185),
                      xi_189),
                  xi_236));
          _mm256_store_pd(
              &_data_pdfs_20_310_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(_mm256_add_pd(forceTerm_10, xi_175),
                                    xi_187),
                      xi_189),
                  xi_225));
          _mm256_store_pd(
              &_data_pdfs_20_311_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(_mm256_add_pd(forceTerm_11, xi_193),
                                    xi_201),
                      xi_204),
                  xi_220));
          _mm256_store_pd(
              &_data_pdfs_20_312_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(_mm256_add_pd(forceTerm_12, xi_201),
                                    xi_205),
                      xi_206),
                  xi_221));
          _mm256_store_pd(
              &_data_pdfs_20_313_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(_mm256_add_pd(forceTerm_13, xi_210),
                                    xi_211),
                      xi_214),
                  xi_227));
          _mm256_store_pd(
              &_data_pdfs_20_314_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(_mm256_add_pd(forceTerm_14, xi_211),
                                    xi_215),
                      xi_216),
                  xi_235));
          _mm256_store_pd(
              &_data_pdfs_20_315_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(_mm256_add_pd(forceTerm_15, xi_204),
                                    xi_206),
                      xi_218),
                  xi_223));
          _mm256_store_pd(
              &_data_pdfs_20_316_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(_mm256_add_pd(forceTerm_16, xi_193),
                                    xi_205),
                      xi_218),
                  xi_224));
          _mm256_store_pd(
              &_data_pdfs_20_317_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(_mm256_add_pd(forceTerm_17, xi_214),
                                    xi_215),
                      xi_219),
                  xi_240));
          _mm256_store_pd(
              &_data_pdfs_20_318_10[ctr_0],
              _mm256_add_pd(
                  _mm256_add_pd(
                      _mm256_add_pd(_mm256_add_pd(forceTerm_18, xi_210),
                                    xi_216),
                      xi_219),
                  xi_229));
        }
        for (int64_t ctr_0 = (int64_t)((_size_force_0) / (4)) * (4);
             ctr_0 < _size_force_0; ctr_0 += 1) {
          const double xi_220 = _data_pdfs_20_311_10[ctr_0];
          const double xi_221 = _data_pdfs_20_312_10[ctr_0];
          const double xi_222 = _data_pdfs_20_35_10[ctr_0];
          const double xi_223 = _data_pdfs_20_315_10[ctr_0];
          const double xi_224 = _data_pdfs_20_316_10[ctr_0];
          const double xi_225 = _data_pdfs_20_310_10[ctr_0];
          const double xi_226 = _data_pdfs_20_37_10[ctr_0];
          const double xi_227 = _data_pdfs_20_313_10[ctr_0];
          const double xi_228 = _data_force_20_31_10[ctr_0];
          const double xi_229 = _data_pdfs_20_318_10[ctr_0];
          const double xi_230 = _data_pdfs_20_30_10[ctr_0];
          const double xi_231 = _data_pdfs_20_31_10[ctr_0];
          const double xi_232 = _data_pdfs_20_32_10[ctr_0];
          const double xi_233 = _data_pdfs_20_38_10[ctr_0];
          const double xi_234 = _data_force_20_32_10[ctr_0];
          const double xi_235 = _data_pdfs_20_314_10[ctr_0];
          const double xi_236 = _data_pdfs_20_39_10[ctr_0];
          const double xi_237 = _data_pdfs_20_36_10[ctr_0];
          const double xi_238 = _data_pdfs_20_33_10[ctr_0];
          const double xi_239 = _data_pdfs_20_34_10[ctr_0];
          const double xi_240 = _data_pdfs_20_317_10[ctr_0];
          const double xi_241 = _data_force_20_30_10[ctr_0];
          const double xi_0 = xi_229 + xi_235;
          const double xi_1 = xi_0 + xi_239;
          const double xi_2 = xi_220 + xi_223 + xi_231;
          const double xi_3 = xi_221 + xi_222;
          const double xi_4 = xi_236 + xi_238;
          const double xi_5 = xi_224 + xi_232;
          const double xi_6 = xi_237 + xi_240;
          const double xi_9 = -xi_236;
          const double xi_10 = -xi_226 + xi_9;
          const double xi_11 = -xi_240;
          const double xi_12 = -xi_227;
          const double xi_13 = -xi_238;
          const double xi_14 = xi_11 + xi_12 + xi_13;
          const double xi_15 = -xi_232;
          const double xi_16 = -xi_225;
          const double xi_17 = xi_15 + xi_16;
          const double xi_18 = -xi_224;
          const double xi_19 = -xi_221;
          const double xi_20 = xi_18 + xi_19;
          const double xi_21 = -xi_229;
          const double xi_22 = xi_11 + xi_21;
          const double xi_23 = -xi_223;
          const double xi_24 = -xi_237;
          const double xi_25 = xi_18 + xi_220 + xi_23 + xi_24;
          const double xi_30 = xi_228 * 0.16666666666666667;
          const double xi_31 = xi_228 * 0.083333333333333333;
          const double xi_43 = xi_241 * 0.16666666666666667;
          const double xi_44 = xi_241 * 0.083333333333333333;
          const double xi_50 = xi_234 * 0.16666666666666667;
          const double xi_51 = xi_234 * 0.083333333333333333;
          const double xi_68 = xi_228 * 0.25000000000000000;
          const double xi_73 = xi_228 * xi_72;
          const double xi_106 = -xi_230;
          const double xi_107 = xi_106 + xi_222 * 3.0 + xi_237 * 3.0;
          const double xi_108 =
              omega_even *
              (xi_107 + xi_220 * -3.0 + xi_221 * -3.0 + xi_223 * -3.0 +
               xi_224 * -3.0 + xi_231 * 3.0 + xi_232 * 3.0);
          const double xi_109 =
              xi_220 * 2.0 + xi_221 * 2.0 + xi_223 * 2.0 + xi_224 * 2.0;
          const double xi_110 = xi_109 + xi_238 * 5.0 + xi_239 * 5.0;
          const double xi_111 =
              omega_even *
              (xi_107 + xi_110 + xi_227 * -5.0 + xi_229 * -5.0 + xi_231 * -2.0 +
               xi_232 * -2.0 + xi_235 * -5.0 + xi_240 * -5.0);
          const double xi_114 = -xi_220;
          const double xi_115 = xi_114 + xi_19;
          const double xi_116 = -xi_233;
          const double xi_119 = -xi_235;
          const double xi_120 = xi_119 + xi_12 + xi_16 + xi_22;
          const double xi_122 = xi_227 * 2.0;
          const double xi_123 = xi_235 * 2.0;
          const double xi_124 = xi_229 * 2.0 + xi_240 * 2.0;
          const double xi_125 =
              omega_even *
              (xi_106 + xi_110 + xi_122 + xi_123 + xi_124 + xi_222 * -4.0 +
               xi_225 * -7.0 + xi_226 * -7.0 + xi_231 * 5.0 + xi_232 * 5.0 +
               xi_233 * -7.0 + xi_236 * -7.0 + xi_237 * -4.0);
          const double xi_126 = xi_114 + xi_221;
          const double xi_127 = xi_126 + xi_15 + xi_224 + xi_23 + xi_231;
          const double xi_129 = xi_127 * xi_128;
          const double xi_130 = xi_226 * 2.0;
          const double xi_131 = xi_225 * 2.0;
          const double xi_132 = xi_233 * -2.0 + xi_236 * 2.0;
          const double xi_133 =
              -xi_130 + xi_131 + xi_132 + xi_15 + xi_2 + xi_20;
          const double xi_135 = xi_133 * xi_134;
          const double xi_136 = -xi_135;
          const double xi_138 = xi_116 + xi_225;
          const double xi_142 = xi_227 + xi_240;
          const double xi_146 = xi_125 * -0.019841269841269841;
          const double xi_154 = xi_119 + xi_227;
          const double xi_155 = xi_13 + xi_154 + xi_21 + xi_239 + xi_240;
          const double xi_156 = xi_128 * xi_155;
          const double xi_157 = xi_1 + xi_130 - xi_131 + xi_132 + xi_14;
          const double xi_158 = xi_134 * xi_157;
          const double xi_160 = -xi_158;
          const double xi_161 = xi_223 + xi_224;
          const double xi_162 = xi_115 + xi_161 + xi_222 + xi_24;
          const double xi_163 = xi_128 * xi_162;
          const double xi_166 = -xi_122 - xi_123 + xi_124 + xi_25 + xi_3;
          const double xi_167 = xi_134 * xi_166;
          const double xi_168 = -xi_167;
          const double xi_170 = xi_167;
          const double xi_174 = xi_125 * 0.013888888888888889;
          const double xi_190 = xi_111 * -0.0071428571428571429;
          const double xi_192 = xi_108 * 0.025000000000000000;
          const double xi_195 = xi_166 * xi_194;
          const double xi_197 = xi_162 * xi_196;
          const double xi_198 = xi_125 * -0.0039682539682539683;
          const double xi_202 = xi_133 * xi_194;
          const double xi_203 = xi_127 * xi_196;
          const double xi_209 = xi_111 * 0.017857142857142857;
          const double xi_212 = xi_155 * xi_196;
          const double xi_213 = xi_157 * xi_194;
          const double xi_32 = rr_0 * xi_31;
          const double xi_45 = rr_0 * xi_44;
          const double xi_52 = rr_0 * xi_51;
          const double xi_55 = xi_241 * xi_54;
          const double xi_60 = xi_228 * xi_54;
          const double xi_82 = xi_234 * xi_54;
          const double vel0Term = xi_1 + xi_225 + xi_233;
          const double vel1Term = xi_2 + xi_226;
          const double vel2Term = xi_227 + xi_3;
          const double rho =
              vel0Term + vel1Term + vel2Term + xi_230 + xi_4 + xi_5 + xi_6;
          const double xi_7 = 1 / (rho);
          const double xi_8 = xi_7 * 0.50000000000000000;
          const double u_0 = xi_241 * xi_8 + xi_7 * (vel0Term + xi_10 + xi_14);
          const double xi_26 = u_0 * xi_241;
          const double xi_38 = xi_26 * 0.16666666666666667;
          const double xi_39 = xi_26 * 0.083333333333333333;
          const double xi_40 = omega_shear * xi_39;
          const double xi_41 = -xi_38 + xi_40;
          const double xi_57 = -xi_26 * xi_56 + xi_38;
          const double xi_58 = -xi_44 + xi_55 + xi_57;
          const double xi_62 = -xi_26 * xi_61;
          const double xi_69 = u_0 * xi_68;
          const double xi_74 = u_0 * xi_73;
          const double xi_78 = xi_44 - xi_55 + xi_57;
          const double xi_85 = -xi_39;
          const double xi_96 = u_0 * xi_234;
          const double xi_97 = xi_96 * 0.25000000000000000;
          const double xi_100 = xi_72 * xi_96;
          const double xi_112 = rho * (u_0 * u_0);
          const double xi_151 = rho * u_0;
          const double xi_152 = -vel0Term + xi_142 + xi_151 + xi_226 + xi_4;
          const double xi_153 = xi_140 * xi_152;
          const double xi_180 = xi_152 * xi_176;
          const double u_1 =
              xi_228 * xi_8 + xi_7 * (vel1Term + xi_17 + xi_20 + xi_233 + xi_9);
          const double xi_27 = u_1 * xi_228;
          const double xi_33 = xi_27 * 0.16666666666666667;
          const double xi_46 = xi_27 * 0.083333333333333333;
          const double xi_47 = omega_shear * xi_46;
          const double xi_48 = -xi_33 + xi_47;
          const double xi_63 = -xi_27 * xi_61;
          const double xi_70 = u_1 * 0.25000000000000000;
          const double xi_71 = xi_241 * xi_70;
          const double xi_75 = u_1 * xi_72;
          const double xi_76 = xi_241 * xi_75;
          const double xi_77 = -xi_69 - xi_71 + xi_74 + xi_76;
          const double xi_79 = xi_69 + xi_71 - xi_74 - xi_76;
          const double xi_87 = xi_234 * xi_70;
          const double xi_89 = xi_234 * xi_75;
          const double xi_94 = -xi_46;
          const double xi_117 = rho * (u_1 * u_1);
          const double xi_118 = xi_10 + xi_116 + xi_117;
          const double xi_137 = rho * u_1;
          const double xi_139 =
              -vel1Term + xi_137 + xi_138 + xi_221 + xi_236 + xi_5;
          const double xi_141 = xi_139 * xi_140;
          const double xi_172 =
              xi_171 * (u_0 * xi_137 + xi_138 + xi_226 + xi_9);
          const double xi_177 = xi_139 * xi_176;
          const double xi_178 = xi_177;
          const double xi_179 = xi_135 + xi_178;
          const double xi_188 = -xi_177;
          const double xi_189 = xi_136 + xi_188;
          const double xi_204 = xi_178 - xi_202 + xi_203;
          const double xi_205 = xi_188 + xi_202 - xi_203;
          const double u_2 =
              xi_234 * xi_8 + xi_7 * (vel2Term + xi_22 + xi_235 + xi_25);
          const double xi_28 = u_2 * xi_234;
          const double xi_34 = xi_28 * 0.16666666666666667;
          const double xi_35 = xi_28 * 0.083333333333333333;
          const double xi_36 = omega_shear * xi_35;
          const double xi_37 = -xi_34 + xi_36;
          const double xi_42 = -omega_shear * xi_33 +
                               xi_27 * 0.33333333333333333 + xi_37 + xi_41;
          const double xi_49 = -omega_shear * xi_38 +
                               xi_26 * 0.33333333333333333 + xi_37 + xi_48;
          const double xi_53 = -omega_shear * xi_34 +
                               xi_28 * 0.33333333333333333 + xi_41 + xi_48;
          const double xi_59 = -xi_35;
          const double xi_64 = -xi_28 * xi_61;
          const double xi_65 = -xi_27 * xi_56 + xi_33 + xi_62 + xi_63 + xi_64;
          const double xi_66 = xi_31 - xi_60 + xi_65;
          const double xi_67 = xi_36 + xi_59 + xi_66;
          const double xi_80 = -xi_31 + xi_60 + xi_65;
          const double xi_81 = xi_36 + xi_59 + xi_80;
          const double xi_83 = -xi_28 * xi_56 + xi_34;
          const double xi_84 = xi_51 - xi_82 + xi_83;
          const double xi_86 = xi_40 + xi_66 + xi_85;
          const double xi_88 = u_2 * xi_68;
          const double xi_90 = u_2 * xi_73;
          const double xi_91 = xi_87 + xi_88 - xi_89 - xi_90;
          const double xi_92 = xi_40 + xi_80 + xi_85;
          const double xi_93 = -xi_87 - xi_88 + xi_89 + xi_90;
          const double xi_95 = xi_47 + xi_62 + xi_63 + xi_64 + xi_84 + xi_94;
          const double xi_98 = u_2 * xi_241;
          const double xi_99 = xi_98 * 0.25000000000000000;
          const double xi_101 = xi_72 * xi_98;
          const double xi_102 = xi_100 + xi_101 - xi_97 - xi_99;
          const double xi_103 = -xi_100 - xi_101 + xi_97 + xi_99;
          const double xi_104 = -xi_51 + xi_82 + xi_83;
          const double xi_105 = xi_104 + xi_47 + xi_62 + xi_63 + xi_64 + xi_94;
          const double xi_113 = rho * (u_2 * u_2);
          const double xi_121 =
              omega_bulk * (xi_112 + xi_113 + xi_115 + xi_118 + xi_120 + xi_18 +
                            xi_23 + xi_230);
          const double xi_143 = -xi_113 + xi_222 + xi_237;
          const double xi_144 =
              omega_shear * (xi_0 + xi_118 + xi_142 + xi_143 + xi_17 - xi_231);
          const double xi_145 = xi_144 * 0.12500000000000000;
          const double xi_147 =
              omega_shear *
              (xi_10 + xi_109 + xi_112 * 2.0 + xi_116 - xi_117 + xi_120 +
               xi_143 + xi_231 + xi_232 + xi_238 * -2.0 + xi_239 * -2.0);
          const double xi_148 = xi_147 * -0.041666666666666667;
          const double xi_149 = xi_108 * -0.050000000000000000 + xi_148;
          const double xi_150 =
              xi_111 * 0.014285714285714286 + xi_145 + xi_146 + xi_149;
          const double xi_159 = xi_111 * -0.035714285714285714 + xi_146 +
                                xi_147 * 0.083333333333333333;
          const double xi_164 =
              rho * u_2 - vel2Term + xi_114 + xi_119 + xi_161 + xi_229 + xi_6;
          const double xi_165 = xi_140 * xi_164;
          const double xi_169 = xi_111 * -0.021428571428571429 +
                                xi_125 * 0.015873015873015873 - xi_145 + xi_149;
          const double xi_173 = xi_144 * 0.062500000000000000;
          const double xi_175 = -xi_172 + xi_173 + xi_174;
          const double xi_181 = xi_121 * 0.041666666666666667;
          const double xi_182 = xi_147 * 0.020833333333333333 + xi_181;
          const double xi_183 = -xi_180 + xi_182;
          const double xi_184 = xi_160 + xi_183;
          const double xi_185 = xi_172 + xi_173 + xi_174;
          const double xi_186 = xi_180 + xi_182;
          const double xi_187 = xi_158 + xi_186;
          const double xi_191 =
              xi_171 * (u_2 * xi_137 + xi_126 + xi_18 + xi_223);
          const double xi_193 = xi_148 + xi_181 + xi_190 + xi_191 + xi_192;
          const double xi_199 = xi_164 * xi_176;
          const double xi_200 = xi_198 + xi_199;
          const double xi_201 = -xi_195 + xi_197 + xi_200;
          const double xi_206 = xi_148 + xi_181 + xi_190 - xi_191 + xi_192;
          const double xi_207 =
              xi_171 * (u_2 * xi_151 + xi_11 + xi_154 + xi_229);
          const double xi_208 = -xi_173;
          const double xi_210 = -xi_207 + xi_208 + xi_209;
          const double xi_211 = xi_170 + xi_200;
          const double xi_214 = xi_183 - xi_212 + xi_213;
          const double xi_215 = xi_207 + xi_208 + xi_209;
          const double xi_216 = xi_186 + xi_212 - xi_213;
          const double xi_217 = xi_198 - xi_199;
          const double xi_218 = xi_195 - xi_197 + xi_217;
          const double xi_219 = xi_168 + xi_217;
          const double forceTerm_0 = xi_26 * xi_29 - xi_26 + xi_27 * xi_29 -
                                     xi_27 + xi_28 * xi_29 - xi_28;
          const double forceTerm_1 = xi_30 - xi_32 + xi_42;
          const double forceTerm_2 = -xi_30 + xi_32 + xi_42;
          const double forceTerm_3 = -xi_43 + xi_45 + xi_49;
          const double forceTerm_4 = xi_43 - xi_45 + xi_49;
          const double forceTerm_5 = xi_50 - xi_52 + xi_53;
          const double forceTerm_6 = -xi_50 + xi_52 + xi_53;
          const double forceTerm_7 = xi_58 + xi_67 + xi_77;
          const double forceTerm_8 = xi_67 + xi_78 + xi_79;
          const double forceTerm_9 = xi_58 + xi_79 + xi_81;
          const double forceTerm_10 = xi_77 + xi_78 + xi_81;
          const double forceTerm_11 = xi_84 + xi_86 + xi_91;
          const double forceTerm_12 = xi_84 + xi_92 + xi_93;
          const double forceTerm_13 = xi_102 + xi_58 + xi_95;
          const double forceTerm_14 = xi_103 + xi_78 + xi_95;
          const double forceTerm_15 = xi_104 + xi_86 + xi_93;
          const double forceTerm_16 = xi_104 + xi_91 + xi_92;
          const double forceTerm_17 = xi_103 + xi_105 + xi_58;
          const double forceTerm_18 = xi_102 + xi_105 + xi_78;
          _data_pdfs_20_30_10[ctr_0] =
              forceTerm_0 + xi_108 * 0.10000000000000000 +
              xi_111 * 0.042857142857142857 + xi_121 * -0.50000000000000000 +
              xi_125 * 0.023809523809523810 + xi_230;
          _data_pdfs_20_31_10[ctr_0] =
              forceTerm_1 - xi_129 + xi_136 + xi_141 + xi_150 + xi_231;
          _data_pdfs_20_32_10[ctr_0] =
              forceTerm_2 + xi_129 + xi_135 - xi_141 + xi_150 + xi_232;
          _data_pdfs_20_33_10[ctr_0] =
              forceTerm_3 - xi_153 + xi_156 + xi_158 + xi_159 + xi_238;
          _data_pdfs_20_34_10[ctr_0] =
              forceTerm_4 + xi_153 - xi_156 + xi_159 + xi_160 + xi_239;
          _data_pdfs_20_35_10[ctr_0] =
              forceTerm_5 - xi_163 + xi_165 + xi_168 + xi_169 + xi_222;
          _data_pdfs_20_36_10[ctr_0] =
              forceTerm_6 + xi_163 - xi_165 + xi_169 + xi_170 + xi_237;
          _data_pdfs_20_37_10[ctr_0] =
              forceTerm_7 + xi_175 + xi_179 + xi_184 + xi_226;
          _data_pdfs_20_38_10[ctr_0] =
              forceTerm_8 + xi_179 + xi_185 + xi_187 + xi_233;
          _data_pdfs_20_39_10[ctr_0] =
              forceTerm_9 + xi_184 + xi_185 + xi_189 + xi_236;
          _data_pdfs_20_310_10[ctr_0] =
              forceTerm_10 + xi_175 + xi_187 + xi_189 + xi_225;
          _data_pdfs_20_311_10[ctr_0] =
              forceTerm_11 + xi_193 + xi_201 + xi_204 + xi_220;
          _data_pdfs_20_312_10[ctr_0] =
              forceTerm_12 + xi_201 + xi_205 + xi_206 + xi_221;
          _data_pdfs_20_313_10[ctr_0] =
              forceTerm_13 + xi_210 + xi_211 + xi_214 + xi_227;
          _data_pdfs_20_314_10[ctr_0] =
              forceTerm_14 + xi_211 + xi_215 + xi_216 + xi_235;
          _data_pdfs_20_315_10[ctr_0] =
              forceTerm_15 + xi_204 + xi_206 + xi_218 + xi_223;
          _data_pdfs_20_316_10[ctr_0] =
              forceTerm_16 + xi_193 + xi_205 + xi_218 + xi_224;
          _data_pdfs_20_317_10[ctr_0] =
              forceTerm_17 + xi_214 + xi_215 + xi_219 + xi_240;
          _data_pdfs_20_318_10[ctr_0] =
              forceTerm_18 + xi_210 + xi_216 + xi_219 + xi_229;
        }
      }
    }
  }
}
} // namespace internal_87ff45d0a5307288e29e018cd62db61e

void CollideSweepDoublePrecisionAVX::run(IBlock *block) {
  auto force = block->getData<field::GhostLayerField<double, 3>>(forceID);
  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);

  auto &omega_bulk = this->omega_bulk_;
  auto &omega_even = this->omega_even_;
  auto &omega_shear = this->omega_shear_;
  auto &omega_odd = this->omega_odd_;
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force->nrOfGhostLayers()));
  double *RESTRICT const _data_force = force->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()));
  double *RESTRICT _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(),
                                int64_t(cell_idx_c(force->xSize()) + 0));
  const int64_t _size_force_0 = int64_t(cell_idx_c(force->xSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(),
                                int64_t(cell_idx_c(force->ySize()) + 0));
  const int64_t _size_force_1 = int64_t(cell_idx_c(force->ySize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(),
                                int64_t(cell_idx_c(force->zSize()) + 0));
  const int64_t _size_force_2 = int64_t(cell_idx_c(force->zSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  internal_87ff45d0a5307288e29e018cd62db61e::
      collidesweepdoubleprecisionavx_collidesweepdoubleprecisionavx(
          _data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2,
          _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_1,
          _stride_pdfs_2, _stride_pdfs_3, omega_bulk, omega_even, omega_odd,
          omega_shear);
}

void CollideSweepDoublePrecisionAVX::runOnCellInterval(
    const shared_ptr<StructuredBlockStorage> &blocks,
    const CellInterval &globalCellInterval, cell_idx_t ghostLayers,
    IBlock *block) {
  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto force = block->getData<field::GhostLayerField<double, 3>>(forceID);
  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);

  auto &omega_bulk = this->omega_bulk_;
  auto &omega_even = this->omega_even_;
  auto &omega_shear = this->omega_shear_;
  auto &omega_odd = this->omega_odd_;
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(force->nrOfGhostLayers()));
  double *RESTRICT const _data_force =
      force->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
  double *RESTRICT _data_pdfs =
      pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.xSize()) + 0));
  const int64_t _size_force_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.ySize()) + 0));
  const int64_t _size_force_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.zSize()) + 0));
  const int64_t _size_force_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  internal_87ff45d0a5307288e29e018cd62db61e::
      collidesweepdoubleprecisionavx_collidesweepdoubleprecisionavx(
          _data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2,
          _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_1,
          _stride_pdfs_2, _stride_pdfs_3, omega_bulk, omega_even, omega_odd,
          omega_shear);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif