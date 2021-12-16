// kernel generated with pystencils v0.4.3+4.g30da657, lbmpy v0.4.3+2.g0e17e61,
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
//! \\file CollideSweepSinglePrecisionLeesEdwardsAVX.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

#include <cmath>

#include "CollideSweepSinglePrecisionLeesEdwardsAVX.h"
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

namespace internal_9a18f2f4073cdcc5365cdfddb752069e {
static FUNC_PREFIX void
collidesweepsingleprecisionleesedwardsavx_collidesweepsingleprecisionleesedwardsavx(
    float *RESTRICT const _data_force, float *RESTRICT _data_pdfs,
    float *RESTRICT const _data_velocity, int64_t const _size_force_0,
    int64_t const _size_force_1, int64_t const _size_force_2,
    int64_t const _stride_force_1, int64_t const _stride_force_2,
    int64_t const _stride_force_3, int64_t const _stride_pdfs_1,
    int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3,
    int64_t const _stride_velocity_1, int64_t const _stride_velocity_2,
    int64_t const _stride_velocity_3, float omega_shear) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1) {
    float *RESTRICT _data_pdfs_20_34 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_318 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_35 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_313 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_31 =
        _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    float *RESTRICT _data_velocity_20_32 =
        _data_velocity + _stride_velocity_2 * ctr_2 + 2 * _stride_velocity_3;
    float *RESTRICT _data_pdfs_20_316 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3;
    float *RESTRICT _data_velocity_20_31 =
        _data_velocity + _stride_velocity_2 * ctr_2 + _stride_velocity_3;
    float *RESTRICT _data_velocity_20_30 =
        _data_velocity + _stride_velocity_2 * ctr_2;
    float *RESTRICT _data_pdfs_20_317 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_36 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_32 =
        _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    float *RESTRICT _data_pdfs_20_312 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_33 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    float *RESTRICT _data_pdfs_20_311 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_315 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_314 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    float *RESTRICT _data_pdfs_20_32 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_31 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1) {
      float *RESTRICT _data_pdfs_20_34_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      float *RESTRICT _data_pdfs_20_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_318;
      float *RESTRICT _data_pdfs_20_35_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_35;
      float *RESTRICT _data_pdfs_20_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_313;
      float *RESTRICT _data_pdfs_20_310_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_310;
      float *RESTRICT _data_force_20_31_10 =
          _stride_force_1 * ctr_1 + _data_force_20_31;
      float *RESTRICT _data_velocity_20_32_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_32;
      float *RESTRICT _data_pdfs_20_316_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_316;
      float *RESTRICT _data_velocity_20_31_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_31;
      float *RESTRICT _data_velocity_20_30_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_30;
      float *RESTRICT _data_pdfs_20_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_317;
      float *RESTRICT _data_pdfs_20_38_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_38;
      float *RESTRICT _data_pdfs_20_36_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_36;
      float *RESTRICT _data_pdfs_20_37_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_37;
      float *RESTRICT _data_force_20_32_10 =
          _stride_force_1 * ctr_1 + _data_force_20_32;
      float *RESTRICT _data_pdfs_20_312_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_312;
      float *RESTRICT _data_pdfs_20_33_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      float *RESTRICT _data_pdfs_20_30_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      float *RESTRICT _data_pdfs_20_311_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_311;
      float *RESTRICT _data_pdfs_20_315_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_315;
      float *RESTRICT _data_pdfs_20_39_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_39;
      float *RESTRICT _data_pdfs_20_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_314;
      float *RESTRICT _data_force_20_30_10 =
          _stride_force_1 * ctr_1 + _data_force_20_30;
      float *RESTRICT _data_pdfs_20_32_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_32;
      float *RESTRICT _data_pdfs_20_31_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_31;
      {
        for (int64_t ctr_0 = 0; ctr_0 < (int64_t)((_size_force_0) / (8)) * (8);
             ctr_0 += 8) {
          const __m256 xi_26 = _mm256_load_ps(&_data_pdfs_20_34_10[ctr_0]);
          const __m256 xi_27 = _mm256_load_ps(&_data_pdfs_20_318_10[ctr_0]);
          const __m256 xi_28 = _mm256_load_ps(&_data_pdfs_20_35_10[ctr_0]);
          const __m256 xi_29 = _mm256_load_ps(&_data_pdfs_20_313_10[ctr_0]);
          const __m256 xi_30 = _mm256_load_ps(&_data_pdfs_20_310_10[ctr_0]);
          const __m256 xi_31 = _mm256_load_ps(&_data_force_20_31_10[ctr_0]);
          const __m256 xi_32 = _mm256_load_ps(&_data_velocity_20_32_10[ctr_0]);
          const __m256 xi_33 = _mm256_load_ps(&_data_pdfs_20_316_10[ctr_0]);
          const __m256 xi_34 = _mm256_load_ps(&_data_velocity_20_31_10[ctr_0]);
          const __m256 xi_35 = _mm256_load_ps(&_data_velocity_20_30_10[ctr_0]);
          const __m256 xi_36 = _mm256_load_ps(&_data_pdfs_20_317_10[ctr_0]);
          const __m256 xi_37 = _mm256_load_ps(&_data_pdfs_20_38_10[ctr_0]);
          const __m256 xi_38 = _mm256_load_ps(&_data_pdfs_20_36_10[ctr_0]);
          const __m256 xi_39 = _mm256_load_ps(&_data_pdfs_20_37_10[ctr_0]);
          const __m256 xi_40 = _mm256_load_ps(&_data_force_20_32_10[ctr_0]);
          const __m256 xi_41 = _mm256_load_ps(&_data_pdfs_20_312_10[ctr_0]);
          const __m256 xi_42 = _mm256_load_ps(&_data_pdfs_20_33_10[ctr_0]);
          const __m256 xi_43 = _mm256_load_ps(&_data_pdfs_20_30_10[ctr_0]);
          const __m256 xi_44 = _mm256_load_ps(&_data_pdfs_20_311_10[ctr_0]);
          const __m256 xi_45 = _mm256_load_ps(&_data_pdfs_20_315_10[ctr_0]);
          const __m256 xi_46 = _mm256_load_ps(&_data_pdfs_20_39_10[ctr_0]);
          const __m256 xi_47 = _mm256_load_ps(&_data_pdfs_20_314_10[ctr_0]);
          const __m256 xi_48 = _mm256_load_ps(&_data_force_20_30_10[ctr_0]);
          const __m256 xi_49 = _mm256_load_ps(&_data_pdfs_20_32_10[ctr_0]);
          const __m256 xi_50 = _mm256_load_ps(&_data_pdfs_20_31_10[ctr_0]);
          const __m256 xi_1 = xi_26;
          const __m256 xi_2 = xi_27;
          const __m256 xi_3 = xi_28;
          const __m256 xi_4 = xi_29;
          const __m256 xi_5 = xi_30;
          const __m256 xi_6 = xi_31;
          const __m256 xi_7 = xi_32;
          const __m256 xi_8 = xi_33;
          const __m256 xi_9 = xi_34;
          const __m256 xi_10 = xi_35;
          const __m256 xi_11 = xi_36;
          const __m256 xi_12 = xi_37;
          const __m256 xi_13 = xi_38;
          const __m256 xi_14 = xi_39;
          const __m256 xi_15 = xi_40;
          const __m256 xi_16 = xi_41;
          const __m256 xi_17 = xi_42;
          const __m256 xi_18 = xi_43;
          const __m256 xi_19 = xi_45;
          const __m256 xi_20 = xi_46;
          const __m256 xi_21 = xi_48;
          const __m256 xi_22 = xi_47;
          const __m256 xi_23 = xi_44;
          const __m256 xi_24 = xi_49;
          const __m256 xi_25 = xi_50;
          const __m256 rho = _mm256_add_ps(
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_add_ps(
                                                      _mm256_add_ps(
                                                          _mm256_add_ps(
                                                              _mm256_add_ps(
                                                                  _mm256_add_ps(
                                                                      _mm256_add_ps(
                                                                          _mm256_add_ps(
                                                                              _mm256_add_ps(
                                                                                  xi_1,
                                                                                  xi_11),
                                                                              xi_12),
                                                                          xi_13),
                                                                      xi_14),
                                                                  xi_16),
                                                              xi_17),
                                                          xi_18),
                                                      xi_19),
                                                  xi_2),
                                              xi_20),
                                          xi_22),
                                      xi_23),
                                  xi_24),
                              xi_25),
                          xi_3),
                      xi_4),
                  xi_5),
              xi_8);
          const __m256 u1Pu2 = _mm256_add_ps(xi_7, xi_9);
          const __m256 u1Mu2 = _mm256_add_ps(
              _mm256_mul_ps(xi_7, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)),
              xi_9);
          _mm256_store_ps(
              &_data_pdfs_20_30_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(_mm256_mul_ps(xi_10, xi_21),
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0)),
                                  _mm256_mul_ps(_mm256_mul_ps(xi_15, xi_7),
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0))),
                              _mm256_mul_ps(_mm256_mul_ps(xi_6, xi_9),
                                            _mm256_set_ps(-1.0, -1.0, -1.0,
                                                          -1.0, -1.0, -1.0,
                                                          -1.0, -1.0))),
                          _mm256_set_ps(omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f)),
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_mul_ps(
                                              rho, _mm256_set_ps(
                                                       0.333333333333333f,
                                                       0.333333333333333f,
                                                       0.333333333333333f,
                                                       0.333333333333333f,
                                                       0.333333333333333f,
                                                       0.333333333333333f,
                                                       0.333333333333333f,
                                                       0.333333333333333f)),
                                          _mm256_mul_ps(
                                              xi_18,
                                              _mm256_set_ps(-1.0, -1.0, -1.0,
                                                            -1.0, -1.0, -1.0,
                                                            -1.0, -1.0))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              rho,
                                              (_mm256_mul_ps(xi_10, xi_10))),
                                          _mm256_set_ps(-0.333333333333333f,
                                                        -0.333333333333333f,
                                                        -0.333333333333333f,
                                                        -0.333333333333333f,
                                                        -0.333333333333333f,
                                                        -0.333333333333333f,
                                                        -0.333333333333333f,
                                                        -0.333333333333333f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          rho, (_mm256_mul_ps(xi_7, xi_7))),
                                      _mm256_set_ps(-0.333333333333333f,
                                                    -0.333333333333333f,
                                                    -0.333333333333333f,
                                                    -0.333333333333333f,
                                                    -0.333333333333333f,
                                                    -0.333333333333333f,
                                                    -0.333333333333333f,
                                                    -0.333333333333333f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(rho,
                                                (_mm256_mul_ps(xi_9, xi_9))),
                                  _mm256_set_ps(
                                      -0.333333333333333f, -0.333333333333333f,
                                      -0.333333333333333f, -0.333333333333333f,
                                      -0.333333333333333f, -0.333333333333333f,
                                      -0.333333333333333f,
                                      -0.333333333333333f))),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  xi_18));
          _mm256_store_ps(
              &_data_pdfs_20_31_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_6,
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  xi_9,
                                                  _mm256_set_ps(
                                                      2.0f, 2.0f, 2.0f, 2.0f,
                                                      2.0f, 2.0f, 2.0f, 2.0f)),
                                              _mm256_set_ps(1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f))),
                                      _mm256_set_ps(0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f)),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_21,
                                          _mm256_add_ps(
                                              xi_10,
                                              _mm256_set_ps(
                                                  ((ctr_1 >= 63) ? (-1.0f)
                                                                 : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 >= 63) ? (-1.0f)
                                                                 : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 >= 63) ? (-1.0f)
                                                                 : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 >= 63) ? (-1.0f)
                                                                 : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 >= 63) ? (-1.0f)
                                                                 : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 >= 63) ? (-1.0f)
                                                                 : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 >= 63) ? (-1.0f)
                                                                 : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 >= 63) ? (-1.0f)
                                                                 : (0.0f)) *
                                                      0.050000000000000003f))),
                                      _mm256_set_ps(-0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(xi_15, xi_7),
                                  _mm256_set_ps(
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f,
                                      -0.166666666666667f))),
                          _mm256_set_ps(omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f)),
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_add_ps(
                                                      _mm256_add_ps(
                                                          _mm256_mul_ps(
                                                              (_mm256_mul_ps(
                                                                  xi_9, xi_9)),
                                                              _mm256_set_ps(
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f)),
                                                          _mm256_mul_ps(
                                                              xi_9,
                                                              _mm256_set_ps(
                                                                  0.166666666666667f,
                                                                  0.166666666666667f,
                                                                  0.166666666666667f,
                                                                  0.166666666666667f,
                                                                  0.166666666666667f,
                                                                  0.166666666666667f,
                                                                  0.166666666666667f,
                                                                  0.166666666666667f))),
                                                      _mm256_set_ps(
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f))),
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_set_ps(
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f))),
                                          _mm256_mul_ps(
                                              xi_25,
                                              _mm256_set_ps(-1.0, -1.0, -1.0,
                                                            -1.0, -1.0, -1.0,
                                                            -1.0, -1.0))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              rho,
                                              (_mm256_mul_ps(
                                                  _mm256_add_ps(
                                                      xi_10,
                                                      _mm256_set_ps(
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f)),
                                                  _mm256_add_ps(
                                                      xi_10,
                                                      _mm256_set_ps(
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f))))),
                                          _mm256_set_ps(-0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          rho, (_mm256_mul_ps(xi_7, xi_7))),
                                      _mm256_set_ps(-0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(rho,
                                                (_mm256_mul_ps(xi_9, xi_9))),
                                  _mm256_set_ps(
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f,
                                      -0.166666666666667f))),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  xi_25));
          _mm256_store_ps(
              &_data_pdfs_20_32_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_6,
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  xi_9,
                                                  _mm256_set_ps(
                                                      2.0f, 2.0f, 2.0f, 2.0f,
                                                      2.0f, 2.0f, 2.0f, 2.0f)),
                                              _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                            -1.0f, -1.0f, -1.0f,
                                                            -1.0f, -1.0f))),
                                      _mm256_set_ps(0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f)),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_21,
                                          _mm256_add_ps(
                                              xi_10,
                                              _mm256_set_ps(
                                                  ((ctr_1 <= 0) ? (1.0f)
                                                                : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 <= 0) ? (1.0f)
                                                                : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 <= 0) ? (1.0f)
                                                                : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 <= 0) ? (1.0f)
                                                                : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 <= 0) ? (1.0f)
                                                                : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 <= 0) ? (1.0f)
                                                                : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 <= 0) ? (1.0f)
                                                                : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 <= 0) ? (1.0f)
                                                                : (0.0f)) *
                                                      0.050000000000000003f))),
                                      _mm256_set_ps(-0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(xi_15, xi_7),
                                  _mm256_set_ps(
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f,
                                      -0.166666666666667f))),
                          _mm256_set_ps(omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f)),
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_add_ps(
                                                      _mm256_add_ps(
                                                          _mm256_mul_ps(
                                                              (_mm256_mul_ps(
                                                                  xi_9, xi_9)),
                                                              _mm256_set_ps(
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f)),
                                                          _mm256_mul_ps(
                                                              xi_9,
                                                              _mm256_set_ps(
                                                                  -0.166666666666667f,
                                                                  -0.166666666666667f,
                                                                  -0.166666666666667f,
                                                                  -0.166666666666667f,
                                                                  -0.166666666666667f,
                                                                  -0.166666666666667f,
                                                                  -0.166666666666667f,
                                                                  -0.166666666666667f))),
                                                      _mm256_set_ps(
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f))),
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_set_ps(
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f))),
                                          _mm256_mul_ps(
                                              xi_24,
                                              _mm256_set_ps(-1.0, -1.0, -1.0,
                                                            -1.0, -1.0, -1.0,
                                                            -1.0, -1.0))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              rho,
                                              (_mm256_mul_ps(
                                                  _mm256_add_ps(
                                                      xi_10,
                                                      _mm256_set_ps(
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f)),
                                                  _mm256_add_ps(
                                                      xi_10,
                                                      _mm256_set_ps(
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f))))),
                                          _mm256_set_ps(-0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          rho, (_mm256_mul_ps(xi_7, xi_7))),
                                      _mm256_set_ps(-0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(rho,
                                                (_mm256_mul_ps(xi_9, xi_9))),
                                  _mm256_set_ps(
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f,
                                      -0.166666666666667f))),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  xi_24));
          _mm256_store_ps(
              &_data_pdfs_20_33_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_21,
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  xi_10,
                                                  _mm256_set_ps(
                                                      2.0f, 2.0f, 2.0f, 2.0f,
                                                      2.0f, 2.0f, 2.0f, 2.0f)),
                                              _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                            -1.0f, -1.0f, -1.0f,
                                                            -1.0f, -1.0f))),
                                      _mm256_set_ps(0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f)),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(xi_15, xi_7),
                                      _mm256_set_ps(-0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(xi_6, xi_9),
                                  _mm256_set_ps(
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f,
                                      -0.166666666666667f))),
                          _mm256_set_ps(omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f)),
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_add_ps(
                                                      _mm256_add_ps(
                                                          _mm256_mul_ps(
                                                              (_mm256_mul_ps(
                                                                  xi_10,
                                                                  xi_10)),
                                                              _mm256_set_ps(
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f)),
                                                          _mm256_mul_ps(
                                                              xi_10,
                                                              _mm256_set_ps(
                                                                  -0.166666666666667f,
                                                                  -0.166666666666667f,
                                                                  -0.166666666666667f,
                                                                  -0.166666666666667f,
                                                                  -0.166666666666667f,
                                                                  -0.166666666666667f,
                                                                  -0.166666666666667f,
                                                                  -0.166666666666667f))),
                                                      _mm256_set_ps(
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f))),
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_set_ps(
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f))),
                                          _mm256_mul_ps(
                                              xi_17,
                                              _mm256_set_ps(-1.0, -1.0, -1.0,
                                                            -1.0, -1.0, -1.0,
                                                            -1.0, -1.0))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              rho,
                                              (_mm256_mul_ps(xi_10, xi_10))),
                                          _mm256_set_ps(-0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          rho, (_mm256_mul_ps(xi_7, xi_7))),
                                      _mm256_set_ps(-0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(rho,
                                                (_mm256_mul_ps(xi_9, xi_9))),
                                  _mm256_set_ps(
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f,
                                      -0.166666666666667f))),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  xi_17));
          _mm256_store_ps(
              &_data_pdfs_20_34_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_21,
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  xi_10,
                                                  _mm256_set_ps(
                                                      2.0f, 2.0f, 2.0f, 2.0f,
                                                      2.0f, 2.0f, 2.0f, 2.0f)),
                                              _mm256_set_ps(1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f))),
                                      _mm256_set_ps(0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f)),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(xi_15, xi_7),
                                      _mm256_set_ps(-0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(xi_6, xi_9),
                                  _mm256_set_ps(
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f,
                                      -0.166666666666667f))),
                          _mm256_set_ps(omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f)),
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_add_ps(
                                                      _mm256_add_ps(
                                                          _mm256_mul_ps(
                                                              (_mm256_mul_ps(
                                                                  xi_10,
                                                                  xi_10)),
                                                              _mm256_set_ps(
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f)),
                                                          _mm256_mul_ps(
                                                              xi_10,
                                                              _mm256_set_ps(
                                                                  0.166666666666667f,
                                                                  0.166666666666667f,
                                                                  0.166666666666667f,
                                                                  0.166666666666667f,
                                                                  0.166666666666667f,
                                                                  0.166666666666667f,
                                                                  0.166666666666667f,
                                                                  0.166666666666667f))),
                                                      _mm256_set_ps(
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f))),
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_set_ps(
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f))),
                                          _mm256_mul_ps(
                                              xi_1,
                                              _mm256_set_ps(-1.0, -1.0, -1.0,
                                                            -1.0, -1.0, -1.0,
                                                            -1.0, -1.0))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              rho,
                                              (_mm256_mul_ps(xi_10, xi_10))),
                                          _mm256_set_ps(-0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          rho, (_mm256_mul_ps(xi_7, xi_7))),
                                      _mm256_set_ps(-0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(rho,
                                                (_mm256_mul_ps(xi_9, xi_9))),
                                  _mm256_set_ps(
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f,
                                      -0.166666666666667f))),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  xi_1));
          _mm256_store_ps(
              &_data_pdfs_20_35_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_15,
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  xi_7,
                                                  _mm256_set_ps(
                                                      2.0f, 2.0f, 2.0f, 2.0f,
                                                      2.0f, 2.0f, 2.0f, 2.0f)),
                                              _mm256_set_ps(1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f))),
                                      _mm256_set_ps(0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f)),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(xi_10, xi_21),
                                      _mm256_set_ps(-0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(xi_6, xi_9),
                                  _mm256_set_ps(
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f,
                                      -0.166666666666667f))),
                          _mm256_set_ps(omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f)),
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_add_ps(
                                                      _mm256_add_ps(
                                                          _mm256_mul_ps(
                                                              (_mm256_mul_ps(
                                                                  xi_7, xi_7)),
                                                              _mm256_set_ps(
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f)),
                                                          _mm256_mul_ps(
                                                              xi_7,
                                                              _mm256_set_ps(
                                                                  0.166666666666667f,
                                                                  0.166666666666667f,
                                                                  0.166666666666667f,
                                                                  0.166666666666667f,
                                                                  0.166666666666667f,
                                                                  0.166666666666667f,
                                                                  0.166666666666667f,
                                                                  0.166666666666667f))),
                                                      _mm256_set_ps(
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f))),
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_set_ps(
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f))),
                                          _mm256_mul_ps(
                                              xi_3,
                                              _mm256_set_ps(-1.0, -1.0, -1.0,
                                                            -1.0, -1.0, -1.0,
                                                            -1.0, -1.0))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              rho,
                                              (_mm256_mul_ps(xi_10, xi_10))),
                                          _mm256_set_ps(-0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          rho, (_mm256_mul_ps(xi_7, xi_7))),
                                      _mm256_set_ps(-0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(rho,
                                                (_mm256_mul_ps(xi_9, xi_9))),
                                  _mm256_set_ps(
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f,
                                      -0.166666666666667f))),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  xi_3));
          _mm256_store_ps(
              &_data_pdfs_20_36_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_15,
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  xi_7,
                                                  _mm256_set_ps(
                                                      2.0f, 2.0f, 2.0f, 2.0f,
                                                      2.0f, 2.0f, 2.0f, 2.0f)),
                                              _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                            -1.0f, -1.0f, -1.0f,
                                                            -1.0f, -1.0f))),
                                      _mm256_set_ps(0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f,
                                                    0.166666666666667f)),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(xi_10, xi_21),
                                      _mm256_set_ps(-0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(xi_6, xi_9),
                                  _mm256_set_ps(
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f,
                                      -0.166666666666667f))),
                          _mm256_set_ps(omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f)),
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_add_ps(
                                                      _mm256_add_ps(
                                                          _mm256_mul_ps(
                                                              (_mm256_mul_ps(
                                                                  xi_7, xi_7)),
                                                              _mm256_set_ps(
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f,
                                                                  0.333333333333333f)),
                                                          _mm256_mul_ps(
                                                              xi_7,
                                                              _mm256_set_ps(
                                                                  -0.166666666666667f,
                                                                  -0.166666666666667f,
                                                                  -0.166666666666667f,
                                                                  -0.166666666666667f,
                                                                  -0.166666666666667f,
                                                                  -0.166666666666667f,
                                                                  -0.166666666666667f,
                                                                  -0.166666666666667f))),
                                                      _mm256_set_ps(
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f,
                                                          -0.111111111111111f))),
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_set_ps(
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f,
                                                      0.166666666666667f))),
                                          _mm256_mul_ps(
                                              xi_13,
                                              _mm256_set_ps(-1.0, -1.0, -1.0,
                                                            -1.0, -1.0, -1.0,
                                                            -1.0, -1.0))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              rho,
                                              (_mm256_mul_ps(xi_10, xi_10))),
                                          _mm256_set_ps(-0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f,
                                                        -0.166666666666667f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          rho, (_mm256_mul_ps(xi_7, xi_7))),
                                      _mm256_set_ps(-0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f,
                                                    -0.166666666666667f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(rho,
                                                (_mm256_mul_ps(xi_9, xi_9))),
                                  _mm256_set_ps(
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f, -0.166666666666667f,
                                      -0.166666666666667f,
                                      -0.166666666666667f))),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  xi_13));
          _mm256_store_ps(&_data_pdfs_20_37_10[ctr_0],
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  _mm256_mul_ps(
                                                      xi_6,
                                                      _mm256_add_ps(
                                                          _mm256_add_ps(
                                                              _mm256_add_ps(
                                                                  _mm256_mul_ps(
                                                                      xi_10,
                                                                      _mm256_set_ps(
                                                                          -3.0f,
                                                                          -3.0f,
                                                                          -3.0f,
                                                                          -3.0f,
                                                                          -3.0f,
                                                                          -3.0f,
                                                                          -3.0f,
                                                                          -3.0f)),
                                                                  _mm256_mul_ps(
                                                                      xi_9,
                                                                      _mm256_set_ps(
                                                                          2.0f,
                                                                          2.0f,
                                                                          2.0f,
                                                                          2.0f,
                                                                          2.0f,
                                                                          2.0f,
                                                                          2.0f,
                                                                          2.0f))),
                                                              _mm256_set_ps(
                                                                  ((ctr_1 >= 63)
                                                                       ? (-1.0f)
                                                                       : (0.0f)) *
                                                                      -0.15000000000000002f,
                                                                  ((ctr_1 >= 63)
                                                                       ? (-1.0f)
                                                                       : (0.0f)) *
                                                                      -0.15000000000000002f,
                                                                  ((ctr_1 >= 63)
                                                                       ? (-1.0f)
                                                                       : (0.0f)) *
                                                                      -0.15000000000000002f,
                                                                  ((ctr_1 >= 63)
                                                                       ? (-1.0f)
                                                                       : (0.0f)) *
                                                                      -0.15000000000000002f,
                                                                  ((ctr_1 >= 63)
                                                                       ? (-1.0f)
                                                                       : (0.0f)) *
                                                                      -0.15000000000000002f,
                                                                  ((ctr_1 >= 63)
                                                                       ? (-1.0f)
                                                                       : (0.0f)) *
                                                                      -0.15000000000000002f,
                                                                  ((ctr_1 >= 63)
                                                                       ? (-1.0f)
                                                                       : (0.0f)) *
                                                                      -0.15000000000000002f,
                                                                  ((ctr_1 >= 63)
                                                                       ? (-1.0f)
                                                                       : (0.0f)) *
                                                                      -0.15000000000000002f)),
                                                          _mm256_set_ps(
                                                              1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f))),
                                                  _mm256_set_ps(
                                                      0.0833333333333333f,
                                                      0.0833333333333333f,
                                                      0.0833333333333333f,
                                                      0.0833333333333333f,
                                                      0.0833333333333333f,
                                                      0.0833333333333333f,
                                                      0.0833333333333333f,
                                                      0.0833333333333333f)),
                                              _mm256_mul_ps(
                                                  _mm256_mul_ps(
                                                      xi_21,
                                                      _mm256_add_ps(
                                                          _mm256_add_ps(
                                                              _mm256_add_ps(
                                                                  _mm256_mul_ps(
                                                                      xi_10,
                                                                      _mm256_set_ps(
                                                                          -2.0f,
                                                                          -2.0f,
                                                                          -2.0f,
                                                                          -2.0f,
                                                                          -2.0f,
                                                                          -2.0f,
                                                                          -2.0f,
                                                                          -2.0f)),
                                                                  _mm256_mul_ps(
                                                                      xi_9,
                                                                      _mm256_set_ps(3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f))),
                                                              _mm256_set_ps(
                                                                  ((ctr_1 >= 63)
                                                                       ? (-1.0f)
                                                                       : (0.0f)) *
                                                                      -0.10000000000000001f,
                                                                  ((ctr_1 >= 63)
                                                                       ? (-1.0f)
                                                                       : (0.0f)) *
                                                                      -0.10000000000000001f,
                                                                  ((ctr_1 >= 63)
                                                                       ? (-1.0f)
                                                                       : (0.0f)) *
                                                                      -0.10000000000000001f,
                                                                  ((ctr_1 >= 63)
                                                                       ? (-1.0f)
                                                                       : (0.0f)) *
                                                                      -0.10000000000000001f,
                                                                  ((ctr_1 >= 63)
                                                                       ? (-1.0f)
                                                                       : (0.0f)) *
                                                                      -0.10000000000000001f,
                                                                  ((ctr_1 >= 63)
                                                                       ? (-1.0f)
                                                                       : (0.0f)) *
                                                                      -0.10000000000000001f,
                                                                  ((ctr_1 >= 63)
                                                                       ? (-1.0f)
                                                                       : (0.0f)) *
                                                                      -0.10000000000000001f,
                                                                  ((ctr_1 >= 63)
                                                                       ? (-1.0f)
                                                                       : (0.0f)) *
                                                                      -0.10000000000000001f)),
                                                          _mm256_set_ps(
                                                              1.0f,
                                                              1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f))),
                                                  _mm256_set_ps(
                                                      -0.0833333333333333f,
                                                      -0.0833333333333333f,
                                                      -0.0833333333333333f,
                                                      -0.0833333333333333f,
                                                      -0.0833333333333333f,
                                                      -0.0833333333333333f,
                                                      -0.0833333333333333f,
                                                      -0.0833333333333333f))),
                                          _mm256_mul_ps(
                                              _mm256_mul_ps(xi_15, xi_7),
                                              _mm256_set_ps(
                                                  -0.0833333333333333f,
                                                  -0.0833333333333333f,
                                                  -0.0833333333333333f,
                                                  -0.0833333333333333f,
                                                  -0.0833333333333333f,
                                                  -0.0833333333333333f,
                                                  -0.0833333333333333f,
                                                  -0.0833333333333333f))),
                                      _mm256_set_ps(omega_shear * -0.5f + 1.0f,
                                                    omega_shear * -0.5f + 1.0f,
                                                    omega_shear * -0.5f + 1.0f,
                                                    omega_shear * -0.5f + 1.0f,
                                                    omega_shear * -0.5f + 1.0f,
                                                    omega_shear * -0.5f + 1.0f,
                                                    omega_shear * -0.5f + 1.0f,
                                                    omega_shear * -0.5f +
                                                        1.0f)),
                                  _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(rho,
                                                                                                                                    _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps((_mm256_mul_ps(
                                                                                                                                                                                                                            _mm256_add_ps(
                                                                                                                                                                                                                                _mm256_add_ps(_mm256_mul_ps(xi_9, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0)), xi_10), _mm256_set_ps(((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                                                                                                                                                                                                                                                                                                                                            0.050000000000000003f,
                                                                                                                                                                                                                                                                                                                                                        ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f,
                                                                                                                                                                                                                                                                                                                                                        ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f,
                                                                                                                                                                                                                                                                                                                                                        ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f,
                                                                                                                                                                                                                                                                                                                                                        ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f)),
                                                                                                                                                                                                                            _mm256_add_ps(_mm256_add_ps(
                                                                                                                                                                                                                                              _mm256_mul_ps(xi_9, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0)),
                                                                                                                                                                                                                                              xi_10),
                                                                                                                                                                                                                                          _mm256_set_ps(
                                                                                                                                                                                                                                              ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f,
                                                                                                                                                                                                                                              ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f,
                                                                                                                                                                                                                                              ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f)))),
                                                                                                                                                                                                                        _mm256_set_ps(0.125f,
                                                                                                                                                                                                                                      0.125f,
                                                                                                                                                                                                                                      0.125f,
                                                                                                                                                                                                                                      0.125f,
                                                                                                                                                                                                                                      0.125f,
                                                                                                                                                                                                                                      0.125f, 0.125f,
                                                                                                                                                                                                                                      0.125f)),
                                                                                                                                                                                                          _mm256_mul_ps(
                                                                                                                                                                                                              (_mm256_mul_ps(
                                                                                                                                                                                                                  xi_7,
                                                                                                                                                                                                                  xi_7)),
                                                                                                                                                                                                              _mm256_set_ps(
                                                                                                                                                                                                                  0.0416666666666667f, 0.0416666666666667f,
                                                                                                                                                                                                                  0.0416666666666667f,
                                                                                                                                                                                                                  0.0416666666666667f,
                                                                                                                                                                                                                  0.0416666666666667f,
                                                                                                                                                                                                                  0.0416666666666667f,
                                                                                                                                                                                                                  0.0416666666666667f,
                                                                                                                                                                                                                  0.0416666666666667f))),
                                                                                                                                                                                            _mm256_mul_ps(
                                                                                                                                                                                                xi_10, _mm256_set_ps(-0.0833333333333333f, -0.0833333333333333f, -0.0833333333333333f,
                                                                                                                                                                                                                     -0.0833333333333333f,
                                                                                                                                                                                                                     -0.0833333333333333f,
                                                                                                                                                                                                                     -0.0833333333333333f,
                                                                                                                                                                                                                     -0.0833333333333333f,
                                                                                                                                                                                                                     -0.0833333333333333f))),
                                                                                                                                                                              _mm256_mul_ps(
                                                                                                                                                                                  xi_9, _mm256_set_ps(0.0833333333333333f, 0.0833333333333333f,
                                                                                                                                                                                                      0.0833333333333333f,
                                                                                                                                                                                                      0.0833333333333333f,
                                                                                                                                                                                                      0.0833333333333333f,
                                                                                                                                                                                                      0.0833333333333333f,
                                                                                                                                                                                                      0.0833333333333333f,
                                                                                                                                                                                                      0.0833333333333333f))),
                                                                                                                                                                _mm256_set_ps(
                                                                                                                                                                    ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * -0.0041666666666666666f,
                                                                                                                                                                    ((ctr_1 >=
                                                                                                                                                                      63)
                                                                                                                                                                         ? (-1.0f)
                                                                                                                                                                         : (0.0f)) *
                                                                                                                                                                        -0.0041666666666666666f,
                                                                                                                                                                    ((ctr_1 >=
                                                                                                                                                                      63)
                                                                                                                                                                         ? (-1.0f)
                                                                                                                                                                         : (0.0f)) *
                                                                                                                                                                        -0.0041666666666666666f,
                                                                                                                                                                    ((
                                                                                                                                                                         ctr_1 >=
                                                                                                                                                                         63)
                                                                                                                                                                         ? (-1.0f)
                                                                                                                                                                         : (0.0f)) *
                                                                                                                                                                        -0.0041666666666666666f,
                                                                                                                                                                    ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * -0.0041666666666666666f, ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * -0.0041666666666666666f, ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * -0.0041666666666666666f,
                                                                                                                                                                    (
                                                                                                                                                                        (ctr_1 >= 63) ? (-1.0f)
                                                                                                                                                                                      : (0.0f)) *
                                                                                                                                                                        -0.0041666666666666666f)),
                                                                                                                                                  _mm256_set_ps(
                                                                                                                                                      -0.0138888888888889f,
                                                                                                                                                      -0.0138888888888889f,
                                                                                                                                                      -0.0138888888888889f,
                                                                                                                                                      -0.0138888888888889f,
                                                                                                                                                      -0.0138888888888889f,
                                                                                                                                                      -0.0138888888888889f,
                                                                                                                                                      -0.0138888888888889f,
                                                                                                                                                      -0.0138888888888889f))),
                                                                                                                      _mm256_mul_ps(
                                                                                                                          rho,
                                                                                                                          _mm256_set_ps(
                                                                                                                              0.0416666666666667f,
                                                                                                                              0.0416666666666667f, 0.0416666666666667f,
                                                                                                                              0.0416666666666667f,
                                                                                                                              0.0416666666666667f,
                                                                                                                              0.0416666666666667f,
                                                                                                                              0.0416666666666667f,
                                                                                                                              0.0416666666666667f))),
                                                                                                        _mm256_mul_ps(
                                                                                                            xi_14,
                                                                                                            _mm256_set_ps(
                                                                                                                -1.0,
                                                                                                                -1.0,
                                                                                                                -1.0,
                                                                                                                -1.0,
                                                                                                                -1.0,
                                                                                                                -1.0,
                                                                                                                -1.0,
                                                                                                                -1.0))),
                                                                                          _mm256_mul_ps(_mm256_mul_ps(rho, (_mm256_mul_ps(_mm256_add_ps(xi_10, _mm256_set_ps((
                                                                                                                                                                                 (ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                                                                                                                                                                 0.050000000000000003f,
                                                                                                                                                                             ((ctr_1 >=
                                                                                                                                                                               63)
                                                                                                                                                                                  ? (-1.0f)
                                                                                                                                                                                  : (0.0f)) *
                                                                                                                                                                                 0.050000000000000003f,
                                                                                                                                                                             (
                                                                                                                                                                                 (ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                                                                                                                                                                 0.050000000000000003f,
                                                                                                                                                                             ((ctr_1 >=
                                                                                                                                                                               63)
                                                                                                                                                                                  ? (-1.0f)
                                                                                                                                                                                  : (0.0f)) *
                                                                                                                                                                                 0.050000000000000003f,
                                                                                                                                                                             (
                                                                                                                                                                                 (ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                                                                                                                                                                 0.050000000000000003f,
                                                                                                                                                                             ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f,
                                                                                                                                                                             ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f)),
                                                                                                                                          _mm256_add_ps(xi_10, _mm256_set_ps(((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f,
                                                                                                                                                                             ((ctr_1 >= 63) ? (-1.0f)
                                                                                                                                                                                            : (0.0f)) *
                                                                                                                                                                                 0.050000000000000003f,
                                                                                                                                                                             (
                                                                                                                                                                                 (ctr_1 >= 63) ? (
                                                                                                                                                                                                     -1.0f)
                                                                                                                                                                                               : (0.0f)) *
                                                                                                                                                                                 0.050000000000000003f))))),
                                                                                                        _mm256_set_ps(
                                                                                                            -0.0416666666666667f,
                                                                                                            -0.0416666666666667f,
                                                                                                            -0.0416666666666667f,
                                                                                                            -0.0416666666666667f,
                                                                                                            -0.0416666666666667f,
                                                                                                            -0.0416666666666667f,
                                                                                                            -0.0416666666666667f,
                                                                                                            -0.0416666666666667f))),
                                                                            _mm256_mul_ps(
                                                                                _mm256_mul_ps(
                                                                                    rho,
                                                                                    (_mm256_mul_ps(
                                                                                        xi_7,
                                                                                        xi_7))),
                                                                                _mm256_set_ps(
                                                                                    -0.0416666666666667f,
                                                                                    -0.0416666666666667f,
                                                                                    -0.0416666666666667f,
                                                                                    -0.0416666666666667f,
                                                                                    -0.0416666666666667f,
                                                                                    -0.0416666666666667f,
                                                                                    -0.0416666666666667f,
                                                                                    -0.0416666666666667f))),
                                                              _mm256_mul_ps(
                                                                  _mm256_mul_ps(
                                                                      rho,
                                                                      (_mm256_mul_ps(
                                                                          xi_9,
                                                                          xi_9))),
                                                                  _mm256_set_ps(
                                                                      -0.0416666666666667f,
                                                                      -0.0416666666666667f,
                                                                      -0.0416666666666667f,
                                                                      -0.0416666666666667f,
                                                                      -0.0416666666666667f,
                                                                      -0.0416666666666667f,
                                                                      -0.0416666666666667f,
                                                                      -0.0416666666666667f))),
                                                _mm256_set_ps(
                                                    omega_shear, omega_shear,
                                                    omega_shear, omega_shear,
                                                    omega_shear, omega_shear,
                                                    omega_shear, omega_shear))),
                              xi_14));
          _mm256_store_ps(
              &_data_pdfs_20_38_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_21,
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_add_ps(
                                                      _mm256_mul_ps(
                                                          xi_10,
                                                          _mm256_set_ps(
                                                              2.0f, 2.0f, 2.0f,
                                                              2.0f, 2.0f, 2.0f,
                                                              2.0f, 2.0f)),
                                                      _mm256_mul_ps(
                                                          xi_9,
                                                          _mm256_set_ps(
                                                              3.0f, 3.0f, 3.0f,
                                                              3.0f, 3.0f, 3.0f,
                                                              3.0f, 3.0f))),
                                                  _mm256_set_ps(
                                                      ((ctr_1 >= 63) ? (-1.0f)
                                                                     : (0.0f)) *
                                                          0.10000000000000001f,
                                                      ((ctr_1 >= 63) ? (-1.0f)
                                                                     : (0.0f)) *
                                                          0.10000000000000001f,
                                                      ((ctr_1 >= 63) ? (-1.0f)
                                                                     : (0.0f)) *
                                                          0.10000000000000001f,
                                                      ((ctr_1 >= 63) ? (-1.0f)
                                                                     : (0.0f)) *
                                                          0.10000000000000001f,
                                                      ((ctr_1 >= 63) ? (-1.0f)
                                                                     : (0.0f)) *
                                                          0.10000000000000001f,
                                                      ((ctr_1 >= 63) ? (-1.0f)
                                                                     : (0.0f)) *
                                                          0.10000000000000001f,
                                                      ((ctr_1 >= 63) ? (-1.0f)
                                                                     : (0.0f)) *
                                                          0.10000000000000001f,
                                                      ((ctr_1 >= 63) ? (-1.0f)
                                                                     : (0.0f)) *
                                                          0.10000000000000001f)),
                                              _mm256_set_ps(1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f))),
                                      _mm256_set_ps(0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f)),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_6,
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_add_ps(
                                                      _mm256_mul_ps(
                                                          xi_10,
                                                          _mm256_set_ps(
                                                              3.0f, 3.0f, 3.0f,
                                                              3.0f, 3.0f, 3.0f,
                                                              3.0f, 3.0f)),
                                                      _mm256_mul_ps(
                                                          xi_9,
                                                          _mm256_set_ps(
                                                              2.0f, 2.0f, 2.0f,
                                                              2.0f, 2.0f, 2.0f,
                                                              2.0f, 2.0f))),
                                                  _mm256_set_ps(
                                                      ((ctr_1 >= 63) ? (-1.0f)
                                                                     : (0.0f)) *
                                                          0.15000000000000002f,
                                                      ((ctr_1 >= 63) ? (-1.0f)
                                                                     : (0.0f)) *
                                                          0.15000000000000002f,
                                                      ((ctr_1 >= 63) ? (-1.0f)
                                                                     : (0.0f)) *
                                                          0.15000000000000002f,
                                                      ((ctr_1 >= 63) ? (-1.0f)
                                                                     : (0.0f)) *
                                                          0.15000000000000002f,
                                                      ((ctr_1 >= 63) ? (-1.0f)
                                                                     : (0.0f)) *
                                                          0.15000000000000002f,
                                                      ((ctr_1 >= 63) ? (-1.0f)
                                                                     : (0.0f)) *
                                                          0.15000000000000002f,
                                                      ((ctr_1 >= 63) ? (-1.0f)
                                                                     : (0.0f)) *
                                                          0.15000000000000002f,
                                                      ((ctr_1 >= 63) ? (-1.0f)
                                                                     : (0.0f)) *
                                                          0.15000000000000002f)),
                                              _mm256_set_ps(1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f))),
                                      _mm256_set_ps(0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(xi_15, xi_7),
                                  _mm256_set_ps(-0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f))),
                          _mm256_set_ps(omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f)),
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_add_ps(
                                                      _mm256_add_ps(
                                                          _mm256_add_ps(
                                                              _mm256_add_ps(
                                                                  _mm256_add_ps(
                                                                      _mm256_mul_ps(
                                                                          (_mm256_mul_ps(
                                                                              _mm256_add_ps(
                                                                                  _mm256_add_ps(
                                                                                      xi_10,
                                                                                      xi_9),
                                                                                  _mm256_set_ps(
                                                                                      ((ctr_1 >=
                                                                                        63)
                                                                                           ? (-1.0f)
                                                                                           : (0.0f)) *
                                                                                          0.050000000000000003f,
                                                                                      ((ctr_1 >=
                                                                                        63)
                                                                                           ? (-1.0f)
                                                                                           : (0.0f)) *
                                                                                          0.050000000000000003f,
                                                                                      ((ctr_1 >=
                                                                                        63)
                                                                                           ? (-1.0f)
                                                                                           : (0.0f)) *
                                                                                          0.050000000000000003f,
                                                                                      ((ctr_1 >=
                                                                                        63)
                                                                                           ? (-1.0f)
                                                                                           : (0.0f)) *
                                                                                          0.050000000000000003f,
                                                                                      ((ctr_1 >=
                                                                                        63)
                                                                                           ? (-1.0f)
                                                                                           : (0.0f)) *
                                                                                          0.050000000000000003f,
                                                                                      ((ctr_1 >=
                                                                                        63)
                                                                                           ? (-1.0f)
                                                                                           : (0.0f)) *
                                                                                          0.050000000000000003f,
                                                                                      ((ctr_1 >=
                                                                                        63)
                                                                                           ? (-1.0f)
                                                                                           : (0.0f)) *
                                                                                          0.050000000000000003f,
                                                                                      ((ctr_1 >=
                                                                                        63)
                                                                                           ? (-1.0f)
                                                                                           : (0.0f)) *
                                                                                          0.050000000000000003f)),
                                                                              _mm256_add_ps(
                                                                                  _mm256_add_ps(xi_10, xi_9), _mm256_set_ps(((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) * 0.050000000000000003f)))),
                                                                          _mm256_set_ps(
                                                                              0.125f,
                                                                              0.125f,
                                                                              0.125f,
                                                                              0.125f,
                                                                              0.125f,
                                                                              0.125f,
                                                                              0.125f,
                                                                              0.125f)),
                                                                      _mm256_mul_ps(
                                                                          (_mm256_mul_ps(
                                                                              xi_7,
                                                                              xi_7)),
                                                                          _mm256_set_ps(
                                                                              0.0416666666666667f,
                                                                              0.0416666666666667f,
                                                                              0.0416666666666667f,
                                                                              0.0416666666666667f,
                                                                              0.0416666666666667f,
                                                                              0.0416666666666667f,
                                                                              0.0416666666666667f,
                                                                              0.0416666666666667f))),
                                                                  _mm256_mul_ps(
                                                                      xi_10,
                                                                      _mm256_set_ps(
                                                                          0.0833333333333333f,
                                                                          0.0833333333333333f,
                                                                          0.0833333333333333f,
                                                                          0.0833333333333333f,
                                                                          0.0833333333333333f,
                                                                          0.0833333333333333f,
                                                                          0.0833333333333333f,
                                                                          0.0833333333333333f))),
                                                              _mm256_mul_ps(
                                                                  xi_9,
                                                                  _mm256_set_ps(
                                                                      0.0833333333333333f,
                                                                      0.0833333333333333f,
                                                                      0.0833333333333333f,
                                                                      0.0833333333333333f,
                                                                      0.0833333333333333f,
                                                                      0.0833333333333333f,
                                                                      0.0833333333333333f,
                                                                      0.0833333333333333f))),
                                                          _mm256_set_ps(
                                                              ((ctr_1 >= 63)
                                                                   ? (-1.0f)
                                                                   : (0.0f)) *
                                                                  0.0041666666666666666f,
                                                              ((ctr_1 >= 63)
                                                                   ? (-1.0f)
                                                                   : (0.0f)) *
                                                                  0.0041666666666666666f,
                                                              ((ctr_1 >= 63)
                                                                   ? (-1.0f)
                                                                   : (0.0f)) *
                                                                  0.0041666666666666666f,
                                                              ((ctr_1 >= 63)
                                                                   ? (-1.0f)
                                                                   : (0.0f)) *
                                                                  0.0041666666666666666f,
                                                              ((ctr_1 >= 63)
                                                                   ? (-1.0f)
                                                                   : (0.0f)) *
                                                                  0.0041666666666666666f,
                                                              ((ctr_1 >= 63)
                                                                   ? (-1.0f)
                                                                   : (0.0f)) *
                                                                  0.0041666666666666666f,
                                                              ((ctr_1 >= 63)
                                                                   ? (-1.0f)
                                                                   : (0.0f)) *
                                                                  0.0041666666666666666f,
                                                              ((ctr_1 >= 63)
                                                                   ? (-1.0f)
                                                                   : (0.0f)) *
                                                                  0.0041666666666666666f)),
                                                      _mm256_set_ps(
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f))),
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_set_ps(
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f))),
                                          _mm256_mul_ps(
                                              xi_12,
                                              _mm256_set_ps(-1.0, -1.0, -1.0,
                                                            -1.0, -1.0, -1.0,
                                                            -1.0, -1.0))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              rho,
                                              (_mm256_mul_ps(
                                                  _mm256_add_ps(
                                                      xi_10,
                                                      _mm256_set_ps(
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f)),
                                                  _mm256_add_ps(
                                                      xi_10,
                                                      _mm256_set_ps(
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f))))),
                                          _mm256_set_ps(-0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          rho, (_mm256_mul_ps(xi_7, xi_7))),
                                      _mm256_set_ps(-0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(rho,
                                                (_mm256_mul_ps(xi_9, xi_9))),
                                  _mm256_set_ps(-0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f))),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  xi_12));
          _mm256_store_ps(&_data_pdfs_20_39_10[ctr_0], _mm256_add_ps(
                                                           _mm256_add_ps(
                                                               _mm256_mul_ps(
                                                                   _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(xi_21, _mm256_add_ps(_mm256_add_ps(
                                                                                                                                                    _mm256_add_ps(_mm256_mul_ps(xi_10, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f)), _mm256_mul_ps(xi_9,
                                                                                                                                                                                                                                                                     _mm256_set_ps(3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f))),
                                                                                                                                                    _mm256_set_ps(
                                                                                                                                                        ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.10000000000000001f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.10000000000000001f,
                                                                                                                                                        ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.10000000000000001f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.10000000000000001f,
                                                                                                                                                        ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.10000000000000001f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.10000000000000001f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.10000000000000001f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.10000000000000001f)),
                                                                                                                                                _mm256_set_ps(-1.0f,
                                                                                                                                                              -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))),
                                                                                                             _mm256_set_ps(
                                                                                                                 0.0833333333333333f, 0.0833333333333333f, 0.0833333333333333f, 0.0833333333333333f, 0.0833333333333333f, 0.0833333333333333f, 0.0833333333333333f, 0.0833333333333333f)),
                                                                                               _mm256_mul_ps(_mm256_mul_ps(xi_6, _mm256_add_ps(
                                                                                                                                     _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_10, _mm256_set_ps(3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f)), _mm256_mul_ps(xi_9, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f))), _mm256_set_ps((
                                                                                                                                                                                                                                                                                                                                                             (ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                                                                                                                                                                                                                                                                                                                                             0.15000000000000002f,
                                                                                                                                                                                                                                                                                                                                                         ((ctr_1 <=
                                                                                                                                                                                                                                                                                                                                                           0)
                                                                                                                                                                                                                                                                                                                                                              ? (1.0f)
                                                                                                                                                                                                                                                                                                                                                              : (0.0f)) *
                                                                                                                                                                                                                                                                                                                                                             0.15000000000000002f,
                                                                                                                                                                                                                                                                                                                                                         ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.15000000000000002f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.15000000000000002f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.15000000000000002f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.15000000000000002f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.15000000000000002f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.15000000000000002f)),
                                                                                                                                     _mm256_set_ps(
                                                                                                                                         -1.0f,
                                                                                                                                         -1.0f,
                                                                                                                                         -1.0f,
                                                                                                                                         -1.0f,
                                                                                                                                         -1.0f,
                                                                                                                                         -1.0f,
                                                                                                                                         -1.0f,
                                                                                                                                         -1.0f))),
                                                                                                             _mm256_set_ps(
                                                                                                                 0.0833333333333333f, 0.0833333333333333f, 0.0833333333333333f, 0.0833333333333333f, 0.0833333333333333f, 0.0833333333333333f, 0.0833333333333333f, 0.0833333333333333f))),
                                                                                 _mm256_mul_ps(_mm256_mul_ps(xi_15, xi_7), _mm256_set_ps(-0.0833333333333333f, -0.0833333333333333f, -0.0833333333333333f, -0.0833333333333333f, -0.0833333333333333f, -0.0833333333333333f, -0.0833333333333333f, -0.0833333333333333f))),
                                                                   _mm256_set_ps(
                                                                       omega_shear *
                                                                               -0.5f +
                                                                           1.0f,
                                                                       omega_shear *
                                                                               -0.5f +
                                                                           1.0f,
                                                                       omega_shear *
                                                                               -0.5f +
                                                                           1.0f,
                                                                       omega_shear *
                                                                               -0.5f +
                                                                           1.0f,
                                                                       omega_shear *
                                                                               -0.5f +
                                                                           1.0f,
                                                                       omega_shear *
                                                                               -0.5f +
                                                                           1.0f,
                                                                       omega_shear *
                                                                               -0.5f +
                                                                           1.0f,
                                                                       omega_shear *
                                                                               -0.5f +
                                                                           1.0f)),
                                                               _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(rho, _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(
                                                                                                                                                                                                                    _mm256_add_ps(_mm256_mul_ps((_mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(xi_10, xi_9), _mm256_set_ps(((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f)), _mm256_add_ps(_mm256_add_ps(xi_10, xi_9), _mm256_set_ps(((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ((ctr_1 <=
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   0)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      ? (1.0f)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      : (0.0f)) *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     0.050000000000000003f,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ((ctr_1 <=
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   0)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      ? (1.0f)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      : (0.0f)) *
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     0.050000000000000003f)))),
                                                                                                                                                                                                                                                _mm256_set_ps(0.125f, 0.125f, 0.125f,
                                                                                                                                                                                                                                                              0.125f,
                                                                                                                                                                                                                                                              0.125f,
                                                                                                                                                                                                                                                              0.125f,
                                                                                                                                                                                                                                                              0.125f,
                                                                                                                                                                                                                                                              0.125f)),
                                                                                                                                                                                                                                  _mm256_mul_ps((_mm256_mul_ps(xi_7, xi_7)), _mm256_set_ps(0.0416666666666667f,
                                                                                                                                                                                                                                                                                           0.0416666666666667f, 0.0416666666666667f, 0.0416666666666667f, 0.0416666666666667f, 0.0416666666666667f, 0.0416666666666667f, 0.0416666666666667f))),
                                                                                                                                                                                                                    _mm256_mul_ps(xi_10, _mm256_set_ps(
                                                                                                                                                                                                                                             -0.0833333333333333f, -0.0833333333333333f, -0.0833333333333333f, -0.0833333333333333f, -0.0833333333333333f, -0.0833333333333333f, -0.0833333333333333f, -0.0833333333333333f))),
                                                                                                                                                                                                                _mm256_mul_ps(xi_9,
                                                                                                                                                                                                                              _mm256_set_ps(-0.0833333333333333f,
                                                                                                                                                                                                                                            -0.0833333333333333f, -0.0833333333333333f, -0.0833333333333333f, -0.0833333333333333f,
                                                                                                                                                                                                                                            -0.0833333333333333f, -0.0833333333333333f, -0.0833333333333333f))),
                                                                                                                                                                                                  _mm256_set_ps(
                                                                                                                                                                                                      (
                                                                                                                                                                                                          (
                                                                                                                                                                                                              ctr_1 <=
                                                                                                                                                                                                              0)
                                                                                                                                                                                                              ? (1.0f)
                                                                                                                                                                                                              : (0.0f)) *
                                                                                                                                                                                                          -0.0041666666666666666f,
                                                                                                                                                                                                      ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * -0.0041666666666666666f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * -0.0041666666666666666f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * -0.0041666666666666666f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * -0.0041666666666666666f,
                                                                                                                                                                                                      ((ctr_1 <=
                                                                                                                                                                                                        0)
                                                                                                                                                                                                           ? (1.0f)
                                                                                                                                                                                                           : (0.0f)) *
                                                                                                                                                                                                          -0.0041666666666666666f,
                                                                                                                                                                                                      ((ctr_1 <=
                                                                                                                                                                                                        0)
                                                                                                                                                                                                           ? (1.0f)
                                                                                                                                                                                                           : (0.0f)) *
                                                                                                                                                                                                          -0.0041666666666666666f,
                                                                                                                                                                                                      (
                                                                                                                                                                                                          (ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                                                                                                                                                                                          -0.0041666666666666666f)),
                                                                                                                                                                                    _mm256_set_ps(
                                                                                                                                                                                        -0.0138888888888889f, -0.0138888888888889f, -0.0138888888888889f, -0.0138888888888889f, -0.0138888888888889f, -0.0138888888888889f, -0.0138888888888889f, -0.0138888888888889f))),
                                                                                                                                                   _mm256_mul_ps(
                                                                                                                                                       rho, _mm256_set_ps(
                                                                                                                                                                0.0416666666666667f, 0.0416666666666667f, 0.0416666666666667f, 0.0416666666666667f, 0.0416666666666667f, 0.0416666666666667f, 0.0416666666666667f,
                                                                                                                                                                0.0416666666666667f))),
                                                                                                                                     _mm256_mul_ps(
                                                                                                                                         xi_20,
                                                                                                                                         _mm256_set_ps(
                                                                                                                                             -1.0,
                                                                                                                                             -1.0,
                                                                                                                                             -1.0,
                                                                                                                                             -1.0,
                                                                                                                                             -1.0,
                                                                                                                                             -1.0,
                                                                                                                                             -1.0,
                                                                                                                                             -1.0))),
                                                                                                                       _mm256_mul_ps(
                                                                                                                           _mm256_mul_ps(
                                                                                                                               rho,
                                                                                                                               (
                                                                                                                                   _mm256_mul_ps(_mm256_add_ps(
                                                                                                                                                     xi_10, _mm256_set_ps(((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                                                                                                                                                              0.050000000000000003f,
                                                                                                                                                                          ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f,
                                                                                                                                                                          (
                                                                                                                                                                              (ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                                                                                                                                                              0.050000000000000003f,
                                                                                                                                                                          ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f,
                                                                                                                                                                          ((ctr_1 <= 0)
                                                                                                                                                                               ? (1.0f)
                                                                                                                                                                               : (0.0f)) *
                                                                                                                                                                              0.050000000000000003f,
                                                                                                                                                                          (
                                                                                                                                                                              (ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                                                                                                                                                              0.050000000000000003f)),
                                                                                                                                                 _mm256_add_ps(xi_10, _mm256_set_ps(((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f,
                                                                                                                                                                                    ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                                                                                                                                                                        0.050000000000000003f))))),
                                                                                                                           _mm256_set_ps(
                                                                                                                               -0.0416666666666667f,
                                                                                                                               -0.0416666666666667f, -0.0416666666666667f, -0.0416666666666667f, -0.0416666666666667f,
                                                                                                                               -0.0416666666666667f,
                                                                                                                               -0.0416666666666667f,
                                                                                                                               -0.0416666666666667f))),
                                                                                                         _mm256_mul_ps(
                                                                                                             _mm256_mul_ps(
                                                                                                                 rho,
                                                                                                                 (_mm256_mul_ps(
                                                                                                                     xi_7,
                                                                                                                     xi_7))),
                                                                                                             _mm256_set_ps(
                                                                                                                 -0.0416666666666667f,
                                                                                                                 -0.0416666666666667f,
                                                                                                                 -0.0416666666666667f,
                                                                                                                 -0.0416666666666667f,
                                                                                                                 -0.0416666666666667f,
                                                                                                                 -0.0416666666666667f,
                                                                                                                 -0.0416666666666667f,
                                                                                                                 -0.0416666666666667f))),
                                                                                           _mm256_mul_ps(_mm256_mul_ps(rho, (_mm256_mul_ps(xi_9, xi_9))), _mm256_set_ps(-0.0416666666666667f,
                                                                                                                                                                        -0.0416666666666667f,
                                                                                                                                                                        -0.0416666666666667f,
                                                                                                                                                                        -0.0416666666666667f,
                                                                                                                                                                        -0.0416666666666667f,
                                                                                                                                                                        -0.0416666666666667f,
                                                                                                                                                                        -0.0416666666666667f,
                                                                                                                                                                        -0.0416666666666667f))),
                                                                             _mm256_set_ps(
                                                                                 omega_shear,
                                                                                 omega_shear,
                                                                                 omega_shear,
                                                                                 omega_shear,
                                                                                 omega_shear,
                                                                                 omega_shear,
                                                                                 omega_shear,
                                                                                 omega_shear))),
                                                           xi_20));
          _mm256_store_ps(
              &_data_pdfs_20_310_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_21,
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_add_ps(
                                                      _mm256_mul_ps(
                                                          xi_10,
                                                          _mm256_set_ps(
                                                              2.0f, 2.0f, 2.0f,
                                                              2.0f, 2.0f, 2.0f,
                                                              2.0f, 2.0f)),
                                                      _mm256_mul_ps(
                                                          xi_9,
                                                          _mm256_set_ps(
                                                              -3.0f, -3.0f,
                                                              -3.0f, -3.0f,
                                                              -3.0f, -3.0f,
                                                              -3.0f, -3.0f))),
                                                  _mm256_set_ps(
                                                      ((ctr_1 <= 0) ? (1.0f)
                                                                    : (0.0f)) *
                                                          0.10000000000000001f,
                                                      ((ctr_1 <= 0) ? (1.0f)
                                                                    : (0.0f)) *
                                                          0.10000000000000001f,
                                                      ((ctr_1 <= 0) ? (1.0f)
                                                                    : (0.0f)) *
                                                          0.10000000000000001f,
                                                      ((ctr_1 <= 0) ? (1.0f)
                                                                    : (0.0f)) *
                                                          0.10000000000000001f,
                                                      ((ctr_1 <= 0) ? (1.0f)
                                                                    : (0.0f)) *
                                                          0.10000000000000001f,
                                                      ((ctr_1 <= 0) ? (1.0f)
                                                                    : (0.0f)) *
                                                          0.10000000000000001f,
                                                      ((ctr_1 <= 0) ? (1.0f)
                                                                    : (0.0f)) *
                                                          0.10000000000000001f,
                                                      ((ctr_1 <= 0) ? (1.0f)
                                                                    : (0.0f)) *
                                                          0.10000000000000001f)),
                                              _mm256_set_ps(1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f))),
                                      _mm256_set_ps(0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f)),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_6,
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_add_ps(
                                                      _mm256_mul_ps(
                                                          xi_10,
                                                          _mm256_set_ps(
                                                              3.0f, 3.0f, 3.0f,
                                                              3.0f, 3.0f, 3.0f,
                                                              3.0f, 3.0f)),
                                                      _mm256_mul_ps(
                                                          xi_9,
                                                          _mm256_set_ps(
                                                              -2.0f, -2.0f,
                                                              -2.0f, -2.0f,
                                                              -2.0f, -2.0f,
                                                              -2.0f, -2.0f))),
                                                  _mm256_set_ps(
                                                      ((ctr_1 <= 0) ? (1.0f)
                                                                    : (0.0f)) *
                                                          0.15000000000000002f,
                                                      ((ctr_1 <= 0) ? (1.0f)
                                                                    : (0.0f)) *
                                                          0.15000000000000002f,
                                                      ((ctr_1 <= 0) ? (1.0f)
                                                                    : (0.0f)) *
                                                          0.15000000000000002f,
                                                      ((ctr_1 <= 0) ? (1.0f)
                                                                    : (0.0f)) *
                                                          0.15000000000000002f,
                                                      ((ctr_1 <= 0) ? (1.0f)
                                                                    : (0.0f)) *
                                                          0.15000000000000002f,
                                                      ((ctr_1 <= 0) ? (1.0f)
                                                                    : (0.0f)) *
                                                          0.15000000000000002f,
                                                      ((ctr_1 <= 0) ? (1.0f)
                                                                    : (0.0f)) *
                                                          0.15000000000000002f,
                                                      ((ctr_1 <= 0) ? (1.0f)
                                                                    : (0.0f)) *
                                                          0.15000000000000002f)),
                                              _mm256_set_ps(1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f))),
                                      _mm256_set_ps(-0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(xi_15, xi_7),
                                  _mm256_set_ps(-0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f))),
                          _mm256_set_ps(omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f)),
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  rho, _mm256_add_ps(
                                                           _mm256_add_ps(
                                                               _mm256_add_ps(
                                                                   _mm256_add_ps(
                                                                       _mm256_add_ps(_mm256_mul_ps((
                                                                                                       _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_9, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0)), xi_10), _mm256_set_ps(((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f,
                                                                                                                                                                                                                                                           (
                                                                                                                                                                                                                                                               (
                                                                                                                                                                                                                                                                   ctr_1 <=
                                                                                                                                                                                                                                                                   0)
                                                                                                                                                                                                                                                                   ? (1.0f)
                                                                                                                                                                                                                                                                   : (0.0f)) *
                                                                                                                                                                                                                                                               0.050000000000000003f,
                                                                                                                                                                                                                                                           ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f,
                                                                                                                                                                                                                                                           (
                                                                                                                                                                                                                                                               (ctr_1 <= 0) ? (1.0f)
                                                                                                                                                                                                                                                                            : (0.0f)) *
                                                                                                                                                                                                                                                               0.050000000000000003f)),
                                                                                                                     _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_9, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0)), xi_10), _mm256_set_ps(((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f, ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.050000000000000003f,
                                                                                                                                                                                                                                                           (
                                                                                                                                                                                                                                                               (
                                                                                                                                                                                                                                                                   ctr_1 <=
                                                                                                                                                                                                                                                                   0)
                                                                                                                                                                                                                                                                   ? (1.0f)
                                                                                                                                                                                                                                                                   : (0.0f)) *
                                                                                                                                                                                                                                                               0.050000000000000003f)))),
                                                                                                   _mm256_set_ps(
                                                                                                       0.125f,
                                                                                                       0.125f,
                                                                                                       0.125f,
                                                                                                       0.125f,
                                                                                                       0.125f,
                                                                                                       0.125f,
                                                                                                       0.125f,
                                                                                                       0.125f)),
                                                                                     _mm256_mul_ps(
                                                                                         (_mm256_mul_ps(
                                                                                             xi_7,
                                                                                             xi_7)),
                                                                                         _mm256_set_ps(
                                                                                             0.0416666666666667f,
                                                                                             0.0416666666666667f,
                                                                                             0.0416666666666667f,
                                                                                             0.0416666666666667f,
                                                                                             0.0416666666666667f,
                                                                                             0.0416666666666667f,
                                                                                             0.0416666666666667f,
                                                                                             0.0416666666666667f))),
                                                                       _mm256_mul_ps(
                                                                           xi_10,
                                                                           _mm256_set_ps(
                                                                               0.0833333333333333f,
                                                                               0.0833333333333333f,
                                                                               0.0833333333333333f,
                                                                               0.0833333333333333f,
                                                                               0.0833333333333333f,
                                                                               0.0833333333333333f,
                                                                               0.0833333333333333f,
                                                                               0.0833333333333333f))),
                                                                   _mm256_mul_ps(
                                                                       xi_9,
                                                                       _mm256_set_ps(
                                                                           -0.0833333333333333f,
                                                                           -0.0833333333333333f,
                                                                           -0.0833333333333333f,
                                                                           -0.0833333333333333f,
                                                                           -0.0833333333333333f,
                                                                           -0.0833333333333333f,
                                                                           -0.0833333333333333f,
                                                                           -0.0833333333333333f))),
                                                               _mm256_set_ps(
                                                                   ((ctr_1 <= 0)
                                                                        ? (1.0f)
                                                                        : (0.0f)) *
                                                                       0.0041666666666666666f,
                                                                   ((ctr_1 <= 0)
                                                                        ? (1.0f)
                                                                        : (0.0f)) *
                                                                       0.0041666666666666666f,
                                                                   ((ctr_1 <= 0)
                                                                        ? (1.0f)
                                                                        : (0.0f)) *
                                                                       0.0041666666666666666f,
                                                                   ((ctr_1 <= 0)
                                                                        ? (1.0f)
                                                                        : (0.0f)) *
                                                                       0.0041666666666666666f,
                                                                   ((ctr_1 <= 0)
                                                                        ? (1.0f)
                                                                        : (0.0f)) *
                                                                       0.0041666666666666666f,
                                                                   ((ctr_1 <= 0)
                                                                        ? (1.0f)
                                                                        : (0.0f)) *
                                                                       0.0041666666666666666f,
                                                                   ((ctr_1 <= 0)
                                                                        ? (1.0f)
                                                                        : (0.0f)) *
                                                                       0.0041666666666666666f,
                                                                   ((ctr_1 <= 0)
                                                                        ? (1.0f)
                                                                        : (0.0f)) *
                                                                       0.0041666666666666666f)),
                                                           _mm256_set_ps(
                                                               -0.0138888888888889f,
                                                               -0.0138888888888889f,
                                                               -0.0138888888888889f,
                                                               -0.0138888888888889f,
                                                               -0.0138888888888889f,
                                                               -0.0138888888888889f,
                                                               -0.0138888888888889f,
                                                               -0.0138888888888889f))),
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_set_ps(
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f))),
                                          _mm256_mul_ps(
                                              xi_5,
                                              _mm256_set_ps(-1.0, -1.0, -1.0,
                                                            -1.0, -1.0, -1.0,
                                                            -1.0, -1.0))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              rho,
                                              (_mm256_mul_ps(
                                                  _mm256_add_ps(
                                                      xi_10,
                                                      _mm256_set_ps(
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f)),
                                                  _mm256_add_ps(
                                                      xi_10,
                                                      _mm256_set_ps(
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f))))),
                                          _mm256_set_ps(-0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          rho, (_mm256_mul_ps(xi_7, xi_7))),
                                      _mm256_set_ps(-0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(rho,
                                                (_mm256_mul_ps(xi_9, xi_9))),
                                  _mm256_set_ps(-0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f))),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  xi_5));
          _mm256_store_ps(
              &_data_pdfs_20_311_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_21,
                                          _mm256_add_ps(
                                              xi_10,
                                              _mm256_set_ps(
                                                  ((ctr_1 >= 63) ? (-1.0f)
                                                                 : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 >= 63) ? (-1.0f)
                                                                 : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 >= 63) ? (-1.0f)
                                                                 : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 >= 63) ? (-1.0f)
                                                                 : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 >= 63) ? (-1.0f)
                                                                 : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 >= 63) ? (-1.0f)
                                                                 : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 >= 63) ? (-1.0f)
                                                                 : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 >= 63) ? (-1.0f)
                                                                 : (0.0f)) *
                                                      0.050000000000000003f))),
                                      _mm256_set_ps(-0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f)),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_15,
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_mul_ps(
                                                      xi_7,
                                                      _mm256_set_ps(
                                                          2.0f, 2.0f, 2.0f,
                                                          2.0f, 2.0f, 2.0f,
                                                          2.0f, 2.0f)),
                                                  _mm256_mul_ps(
                                                      xi_9,
                                                      _mm256_set_ps(
                                                          3.0f, 3.0f, 3.0f,
                                                          3.0f, 3.0f, 3.0f,
                                                          3.0f, 3.0f))),
                                              _mm256_set_ps(1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f))),
                                      _mm256_set_ps(0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(
                                      xi_6,
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  xi_7,
                                                  _mm256_set_ps(
                                                      3.0f, 3.0f, 3.0f, 3.0f,
                                                      3.0f, 3.0f, 3.0f, 3.0f)),
                                              _mm256_mul_ps(
                                                  xi_9,
                                                  _mm256_set_ps(
                                                      2.0f, 2.0f, 2.0f, 2.0f,
                                                      2.0f, 2.0f, 2.0f, 2.0f))),
                                          _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f,
                                                        1.0f, 1.0f, 1.0f,
                                                        1.0f))),
                                  _mm256_set_ps(
                                      0.0833333333333333f, 0.0833333333333333f,
                                      0.0833333333333333f, 0.0833333333333333f,
                                      0.0833333333333333f, 0.0833333333333333f,
                                      0.0833333333333333f,
                                      0.0833333333333333f))),
                          _mm256_set_ps(omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f)),
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_add_ps(
                                                      _mm256_add_ps(
                                                          _mm256_add_ps(
                                                              _mm256_mul_ps(
                                                                  (_mm256_mul_ps(
                                                                      _mm256_add_ps(
                                                                          xi_10,
                                                                          _mm256_set_ps(
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f)),
                                                                      _mm256_add_ps(
                                                                          xi_10,
                                                                          _mm256_set_ps(
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f)))),
                                                                  _mm256_set_ps(
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f)),
                                                              _mm256_mul_ps(
                                                                  (_mm256_mul_ps(
                                                                      u1Pu2,
                                                                      u1Pu2)),
                                                                  _mm256_set_ps(
                                                                      0.125f,
                                                                      0.125f,
                                                                      0.125f,
                                                                      0.125f,
                                                                      0.125f,
                                                                      0.125f,
                                                                      0.125f,
                                                                      0.125f))),
                                                          _mm256_mul_ps(
                                                              u1Pu2,
                                                              _mm256_set_ps(
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f))),
                                                      _mm256_set_ps(
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f))),
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_set_ps(
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f))),
                                          _mm256_mul_ps(
                                              xi_23,
                                              _mm256_set_ps(-1.0, -1.0, -1.0,
                                                            -1.0, -1.0, -1.0,
                                                            -1.0, -1.0))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              rho,
                                              (_mm256_mul_ps(
                                                  _mm256_add_ps(
                                                      xi_10,
                                                      _mm256_set_ps(
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f)),
                                                  _mm256_add_ps(
                                                      xi_10,
                                                      _mm256_set_ps(
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f))))),
                                          _mm256_set_ps(-0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          rho, (_mm256_mul_ps(xi_7, xi_7))),
                                      _mm256_set_ps(-0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(rho,
                                                (_mm256_mul_ps(xi_9, xi_9))),
                                  _mm256_set_ps(-0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f))),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  xi_23));
          _mm256_store_ps(
              &_data_pdfs_20_312_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_21,
                                          _mm256_add_ps(
                                              xi_10,
                                              _mm256_set_ps(
                                                  ((ctr_1 <= 0) ? (1.0f)
                                                                : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 <= 0) ? (1.0f)
                                                                : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 <= 0) ? (1.0f)
                                                                : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 <= 0) ? (1.0f)
                                                                : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 <= 0) ? (1.0f)
                                                                : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 <= 0) ? (1.0f)
                                                                : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 <= 0) ? (1.0f)
                                                                : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 <= 0) ? (1.0f)
                                                                : (0.0f)) *
                                                      0.050000000000000003f))),
                                      _mm256_set_ps(-0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f)),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_15,
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_mul_ps(
                                                      xi_7,
                                                      _mm256_set_ps(
                                                          2.0f, 2.0f, 2.0f,
                                                          2.0f, 2.0f, 2.0f,
                                                          2.0f, 2.0f)),
                                                  _mm256_mul_ps(
                                                      xi_9,
                                                      _mm256_set_ps(
                                                          -3.0f, -3.0f, -3.0f,
                                                          -3.0f, -3.0f, -3.0f,
                                                          -3.0f, -3.0f))),
                                              _mm256_set_ps(1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f))),
                                      _mm256_set_ps(0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(
                                      xi_6,
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  xi_7,
                                                  _mm256_set_ps(
                                                      3.0f, 3.0f, 3.0f, 3.0f,
                                                      3.0f, 3.0f, 3.0f, 3.0f)),
                                              _mm256_mul_ps(
                                                  xi_9, _mm256_set_ps(
                                                            -2.0f, -2.0f, -2.0f,
                                                            -2.0f, -2.0f, -2.0f,
                                                            -2.0f, -2.0f))),
                                          _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f,
                                                        1.0f, 1.0f, 1.0f,
                                                        1.0f))),
                                  _mm256_set_ps(-0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f))),
                          _mm256_set_ps(omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f)),
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_add_ps(
                                                      _mm256_add_ps(
                                                          _mm256_add_ps(
                                                              _mm256_mul_ps(
                                                                  (_mm256_mul_ps(
                                                                      _mm256_add_ps(
                                                                          xi_10,
                                                                          _mm256_set_ps(
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f)),
                                                                      _mm256_add_ps(
                                                                          xi_10,
                                                                          _mm256_set_ps(
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f)))),
                                                                  _mm256_set_ps(
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f)),
                                                              _mm256_mul_ps(
                                                                  (_mm256_mul_ps(
                                                                      u1Mu2,
                                                                      u1Mu2)),
                                                                  _mm256_set_ps(
                                                                      0.125f,
                                                                      0.125f,
                                                                      0.125f,
                                                                      0.125f,
                                                                      0.125f,
                                                                      0.125f,
                                                                      0.125f,
                                                                      0.125f))),
                                                          _mm256_mul_ps(
                                                              u1Mu2,
                                                              _mm256_set_ps(
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f))),
                                                      _mm256_set_ps(
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f))),
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_set_ps(
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f))),
                                          _mm256_mul_ps(
                                              xi_16,
                                              _mm256_set_ps(-1.0, -1.0, -1.0,
                                                            -1.0, -1.0, -1.0,
                                                            -1.0, -1.0))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              rho,
                                              (_mm256_mul_ps(
                                                  _mm256_add_ps(
                                                      xi_10,
                                                      _mm256_set_ps(
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f)),
                                                  _mm256_add_ps(
                                                      xi_10,
                                                      _mm256_set_ps(
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f))))),
                                          _mm256_set_ps(-0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          rho, (_mm256_mul_ps(xi_7, xi_7))),
                                      _mm256_set_ps(-0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(rho,
                                                (_mm256_mul_ps(xi_9, xi_9))),
                                  _mm256_set_ps(-0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f))),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  xi_16));
          _mm256_store_ps(
              &_data_pdfs_20_313_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_15,
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_mul_ps(
                                                      xi_10,
                                                      _mm256_set_ps(
                                                          -3.0f, -3.0f, -3.0f,
                                                          -3.0f, -3.0f, -3.0f,
                                                          -3.0f, -3.0f)),
                                                  _mm256_mul_ps(
                                                      xi_7,
                                                      _mm256_set_ps(
                                                          2.0f, 2.0f, 2.0f,
                                                          2.0f, 2.0f, 2.0f,
                                                          2.0f, 2.0f))),
                                              _mm256_set_ps(1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f))),
                                      _mm256_set_ps(0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f)),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_21,
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_mul_ps(
                                                      xi_10,
                                                      _mm256_set_ps(
                                                          -2.0f, -2.0f, -2.0f,
                                                          -2.0f, -2.0f, -2.0f,
                                                          -2.0f, -2.0f)),
                                                  _mm256_mul_ps(
                                                      xi_7,
                                                      _mm256_set_ps(
                                                          3.0f, 3.0f, 3.0f,
                                                          3.0f, 3.0f, 3.0f,
                                                          3.0f, 3.0f))),
                                              _mm256_set_ps(1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f))),
                                      _mm256_set_ps(-0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(xi_6, xi_9),
                                  _mm256_set_ps(-0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f))),
                          _mm256_set_ps(omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f)),
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_add_ps(
                                                      _mm256_add_ps(
                                                          _mm256_add_ps(
                                                              _mm256_add_ps(
                                                                  _mm256_mul_ps(
                                                                      (_mm256_mul_ps(
                                                                          _mm256_add_ps(
                                                                              _mm256_mul_ps(
                                                                                  xi_7,
                                                                                  _mm256_set_ps(
                                                                                      -1.0,
                                                                                      -1.0,
                                                                                      -1.0,
                                                                                      -1.0,
                                                                                      -1.0,
                                                                                      -1.0,
                                                                                      -1.0,
                                                                                      -1.0)),
                                                                              xi_10),
                                                                          _mm256_add_ps(
                                                                              _mm256_mul_ps(
                                                                                  xi_7,
                                                                                  _mm256_set_ps(
                                                                                      -1.0,
                                                                                      -1.0,
                                                                                      -1.0,
                                                                                      -1.0,
                                                                                      -1.0,
                                                                                      -1.0,
                                                                                      -1.0,
                                                                                      -1.0)),
                                                                              xi_10))),
                                                                      _mm256_set_ps(
                                                                          0.125f,
                                                                          0.125f,
                                                                          0.125f,
                                                                          0.125f,
                                                                          0.125f,
                                                                          0.125f,
                                                                          0.125f,
                                                                          0.125f)),
                                                                  _mm256_mul_ps(
                                                                      (_mm256_mul_ps(
                                                                          xi_9,
                                                                          xi_9)),
                                                                      _mm256_set_ps(
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f))),
                                                              _mm256_mul_ps(
                                                                  xi_10,
                                                                  _mm256_set_ps(
                                                                      -0.0833333333333333f,
                                                                      -0.0833333333333333f,
                                                                      -0.0833333333333333f,
                                                                      -0.0833333333333333f,
                                                                      -0.0833333333333333f,
                                                                      -0.0833333333333333f,
                                                                      -0.0833333333333333f,
                                                                      -0.0833333333333333f))),
                                                          _mm256_mul_ps(
                                                              xi_7,
                                                              _mm256_set_ps(
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f))),
                                                      _mm256_set_ps(
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f))),
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_set_ps(
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f))),
                                          _mm256_mul_ps(
                                              xi_4,
                                              _mm256_set_ps(-1.0, -1.0, -1.0,
                                                            -1.0, -1.0, -1.0,
                                                            -1.0, -1.0))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              rho,
                                              (_mm256_mul_ps(xi_10, xi_10))),
                                          _mm256_set_ps(-0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          rho, (_mm256_mul_ps(xi_7, xi_7))),
                                      _mm256_set_ps(-0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(rho,
                                                (_mm256_mul_ps(xi_9, xi_9))),
                                  _mm256_set_ps(-0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f))),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  xi_4));
          _mm256_store_ps(
              &_data_pdfs_20_314_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_21,
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_mul_ps(
                                                      xi_10,
                                                      _mm256_set_ps(
                                                          2.0f, 2.0f, 2.0f,
                                                          2.0f, 2.0f, 2.0f,
                                                          2.0f, 2.0f)),
                                                  _mm256_mul_ps(
                                                      xi_7,
                                                      _mm256_set_ps(
                                                          3.0f, 3.0f, 3.0f,
                                                          3.0f, 3.0f, 3.0f,
                                                          3.0f, 3.0f))),
                                              _mm256_set_ps(1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f))),
                                      _mm256_set_ps(0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f)),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_15,
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_mul_ps(
                                                      xi_10,
                                                      _mm256_set_ps(
                                                          3.0f, 3.0f, 3.0f,
                                                          3.0f, 3.0f, 3.0f,
                                                          3.0f, 3.0f)),
                                                  _mm256_mul_ps(
                                                      xi_7,
                                                      _mm256_set_ps(
                                                          2.0f, 2.0f, 2.0f,
                                                          2.0f, 2.0f, 2.0f,
                                                          2.0f, 2.0f))),
                                              _mm256_set_ps(1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f))),
                                      _mm256_set_ps(0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(xi_6, xi_9),
                                  _mm256_set_ps(-0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f))),
                          _mm256_set_ps(omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f)),
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_add_ps(
                                                      _mm256_add_ps(
                                                          _mm256_add_ps(
                                                              _mm256_add_ps(
                                                                  _mm256_mul_ps(
                                                                      (_mm256_mul_ps(
                                                                          _mm256_add_ps(
                                                                              xi_10,
                                                                              xi_7),
                                                                          _mm256_add_ps(
                                                                              xi_10,
                                                                              xi_7))),
                                                                      _mm256_set_ps(
                                                                          0.125f,
                                                                          0.125f,
                                                                          0.125f,
                                                                          0.125f,
                                                                          0.125f,
                                                                          0.125f,
                                                                          0.125f,
                                                                          0.125f)),
                                                                  _mm256_mul_ps(
                                                                      (_mm256_mul_ps(
                                                                          xi_9,
                                                                          xi_9)),
                                                                      _mm256_set_ps(
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f))),
                                                              _mm256_mul_ps(
                                                                  xi_10,
                                                                  _mm256_set_ps(
                                                                      0.0833333333333333f,
                                                                      0.0833333333333333f,
                                                                      0.0833333333333333f,
                                                                      0.0833333333333333f,
                                                                      0.0833333333333333f,
                                                                      0.0833333333333333f,
                                                                      0.0833333333333333f,
                                                                      0.0833333333333333f))),
                                                          _mm256_mul_ps(
                                                              xi_7,
                                                              _mm256_set_ps(
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f))),
                                                      _mm256_set_ps(
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f))),
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_set_ps(
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f))),
                                          _mm256_mul_ps(
                                              xi_22,
                                              _mm256_set_ps(-1.0, -1.0, -1.0,
                                                            -1.0, -1.0, -1.0,
                                                            -1.0, -1.0))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              rho,
                                              (_mm256_mul_ps(xi_10, xi_10))),
                                          _mm256_set_ps(-0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          rho, (_mm256_mul_ps(xi_7, xi_7))),
                                      _mm256_set_ps(-0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(rho,
                                                (_mm256_mul_ps(xi_9, xi_9))),
                                  _mm256_set_ps(-0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f))),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  xi_22));
          _mm256_store_ps(
              &_data_pdfs_20_315_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_21,
                                          _mm256_add_ps(
                                              xi_10,
                                              _mm256_set_ps(
                                                  ((ctr_1 >= 63) ? (-1.0f)
                                                                 : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 >= 63) ? (-1.0f)
                                                                 : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 >= 63) ? (-1.0f)
                                                                 : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 >= 63) ? (-1.0f)
                                                                 : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 >= 63) ? (-1.0f)
                                                                 : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 >= 63) ? (-1.0f)
                                                                 : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 >= 63) ? (-1.0f)
                                                                 : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 >= 63) ? (-1.0f)
                                                                 : (0.0f)) *
                                                      0.050000000000000003f))),
                                      _mm256_set_ps(-0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f)),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_6,
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_mul_ps(
                                                      xi_7,
                                                      _mm256_set_ps(
                                                          -3.0f, -3.0f, -3.0f,
                                                          -3.0f, -3.0f, -3.0f,
                                                          -3.0f, -3.0f)),
                                                  _mm256_mul_ps(
                                                      xi_9,
                                                      _mm256_set_ps(
                                                          2.0f, 2.0f, 2.0f,
                                                          2.0f, 2.0f, 2.0f,
                                                          2.0f, 2.0f))),
                                              _mm256_set_ps(1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f))),
                                      _mm256_set_ps(0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(
                                      xi_15,
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  xi_7, _mm256_set_ps(
                                                            -2.0f, -2.0f, -2.0f,
                                                            -2.0f, -2.0f, -2.0f,
                                                            -2.0f, -2.0f)),
                                              _mm256_mul_ps(
                                                  xi_9,
                                                  _mm256_set_ps(
                                                      3.0f, 3.0f, 3.0f, 3.0f,
                                                      3.0f, 3.0f, 3.0f, 3.0f))),
                                          _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f,
                                                        1.0f, 1.0f, 1.0f,
                                                        1.0f))),
                                  _mm256_set_ps(-0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f))),
                          _mm256_set_ps(omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f)),
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_add_ps(
                                                      _mm256_add_ps(
                                                          _mm256_add_ps(
                                                              _mm256_mul_ps(
                                                                  (_mm256_mul_ps(
                                                                      _mm256_add_ps(
                                                                          xi_10,
                                                                          _mm256_set_ps(
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f)),
                                                                      _mm256_add_ps(
                                                                          xi_10,
                                                                          _mm256_set_ps(
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 >=
                                                                                63)
                                                                                   ? (-1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f)))),
                                                                  _mm256_set_ps(
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f)),
                                                              _mm256_mul_ps(
                                                                  (_mm256_mul_ps(
                                                                      u1Mu2,
                                                                      u1Mu2)),
                                                                  _mm256_set_ps(
                                                                      0.125f,
                                                                      0.125f,
                                                                      0.125f,
                                                                      0.125f,
                                                                      0.125f,
                                                                      0.125f,
                                                                      0.125f,
                                                                      0.125f))),
                                                          _mm256_mul_ps(
                                                              u1Mu2,
                                                              _mm256_set_ps(
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f,
                                                                  0.0833333333333333f))),
                                                      _mm256_set_ps(
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f))),
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_set_ps(
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f))),
                                          _mm256_mul_ps(
                                              xi_19,
                                              _mm256_set_ps(-1.0, -1.0, -1.0,
                                                            -1.0, -1.0, -1.0,
                                                            -1.0, -1.0))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              rho,
                                              (_mm256_mul_ps(
                                                  _mm256_add_ps(
                                                      xi_10,
                                                      _mm256_set_ps(
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f)),
                                                  _mm256_add_ps(
                                                      xi_10,
                                                      _mm256_set_ps(
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 >= 63)
                                                               ? (-1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f))))),
                                          _mm256_set_ps(-0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          rho, (_mm256_mul_ps(xi_7, xi_7))),
                                      _mm256_set_ps(-0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(rho,
                                                (_mm256_mul_ps(xi_9, xi_9))),
                                  _mm256_set_ps(-0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f))),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  xi_19));
          _mm256_store_ps(
              &_data_pdfs_20_316_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_21,
                                          _mm256_add_ps(
                                              xi_10,
                                              _mm256_set_ps(
                                                  ((ctr_1 <= 0) ? (1.0f)
                                                                : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 <= 0) ? (1.0f)
                                                                : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 <= 0) ? (1.0f)
                                                                : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 <= 0) ? (1.0f)
                                                                : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 <= 0) ? (1.0f)
                                                                : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 <= 0) ? (1.0f)
                                                                : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 <= 0) ? (1.0f)
                                                                : (0.0f)) *
                                                      0.050000000000000003f,
                                                  ((ctr_1 <= 0) ? (1.0f)
                                                                : (0.0f)) *
                                                      0.050000000000000003f))),
                                      _mm256_set_ps(-0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f)),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_15,
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_mul_ps(
                                                      xi_7,
                                                      _mm256_set_ps(
                                                          2.0f, 2.0f, 2.0f,
                                                          2.0f, 2.0f, 2.0f,
                                                          2.0f, 2.0f)),
                                                  _mm256_mul_ps(
                                                      xi_9,
                                                      _mm256_set_ps(
                                                          3.0f, 3.0f, 3.0f,
                                                          3.0f, 3.0f, 3.0f,
                                                          3.0f, 3.0f))),
                                              _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                            -1.0f, -1.0f, -1.0f,
                                                            -1.0f, -1.0f))),
                                      _mm256_set_ps(0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(
                                      xi_6,
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  xi_7,
                                                  _mm256_set_ps(
                                                      3.0f, 3.0f, 3.0f, 3.0f,
                                                      3.0f, 3.0f, 3.0f, 3.0f)),
                                              _mm256_mul_ps(
                                                  xi_9,
                                                  _mm256_set_ps(
                                                      2.0f, 2.0f, 2.0f, 2.0f,
                                                      2.0f, 2.0f, 2.0f, 2.0f))),
                                          _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                        -1.0f, -1.0f, -1.0f,
                                                        -1.0f, -1.0f))),
                                  _mm256_set_ps(
                                      0.0833333333333333f, 0.0833333333333333f,
                                      0.0833333333333333f, 0.0833333333333333f,
                                      0.0833333333333333f, 0.0833333333333333f,
                                      0.0833333333333333f,
                                      0.0833333333333333f))),
                          _mm256_set_ps(omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f)),
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_add_ps(
                                                      _mm256_add_ps(
                                                          _mm256_add_ps(
                                                              _mm256_mul_ps(
                                                                  (_mm256_mul_ps(
                                                                      _mm256_add_ps(
                                                                          xi_10,
                                                                          _mm256_set_ps(
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f)),
                                                                      _mm256_add_ps(
                                                                          xi_10,
                                                                          _mm256_set_ps(
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f,
                                                                              ((ctr_1 <=
                                                                                0)
                                                                                   ? (1.0f)
                                                                                   : (0.0f)) *
                                                                                  0.050000000000000003f)))),
                                                                  _mm256_set_ps(
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f,
                                                                      0.0416666666666667f)),
                                                              _mm256_mul_ps(
                                                                  (_mm256_mul_ps(
                                                                      u1Pu2,
                                                                      u1Pu2)),
                                                                  _mm256_set_ps(
                                                                      0.125f,
                                                                      0.125f,
                                                                      0.125f,
                                                                      0.125f,
                                                                      0.125f,
                                                                      0.125f,
                                                                      0.125f,
                                                                      0.125f))),
                                                          _mm256_mul_ps(
                                                              u1Pu2,
                                                              _mm256_set_ps(
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f))),
                                                      _mm256_set_ps(
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f))),
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_set_ps(
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f))),
                                          _mm256_mul_ps(
                                              xi_8,
                                              _mm256_set_ps(-1.0, -1.0, -1.0,
                                                            -1.0, -1.0, -1.0,
                                                            -1.0, -1.0))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              rho,
                                              (_mm256_mul_ps(
                                                  _mm256_add_ps(
                                                      xi_10,
                                                      _mm256_set_ps(
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f)),
                                                  _mm256_add_ps(
                                                      xi_10,
                                                      _mm256_set_ps(
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f,
                                                          ((ctr_1 <= 0)
                                                               ? (1.0f)
                                                               : (0.0f)) *
                                                              0.050000000000000003f))))),
                                          _mm256_set_ps(-0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          rho, (_mm256_mul_ps(xi_7, xi_7))),
                                      _mm256_set_ps(-0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(rho,
                                                (_mm256_mul_ps(xi_9, xi_9))),
                                  _mm256_set_ps(-0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f))),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  xi_8));
          _mm256_store_ps(
              &_data_pdfs_20_317_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_21,
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_mul_ps(
                                                      xi_10,
                                                      _mm256_set_ps(
                                                          2.0f, 2.0f, 2.0f,
                                                          2.0f, 2.0f, 2.0f,
                                                          2.0f, 2.0f)),
                                                  _mm256_mul_ps(
                                                      xi_7,
                                                      _mm256_set_ps(
                                                          3.0f, 3.0f, 3.0f,
                                                          3.0f, 3.0f, 3.0f,
                                                          3.0f, 3.0f))),
                                              _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                            -1.0f, -1.0f, -1.0f,
                                                            -1.0f, -1.0f))),
                                      _mm256_set_ps(0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f)),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_15,
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_mul_ps(
                                                      xi_10,
                                                      _mm256_set_ps(
                                                          3.0f, 3.0f, 3.0f,
                                                          3.0f, 3.0f, 3.0f,
                                                          3.0f, 3.0f)),
                                                  _mm256_mul_ps(
                                                      xi_7,
                                                      _mm256_set_ps(
                                                          2.0f, 2.0f, 2.0f,
                                                          2.0f, 2.0f, 2.0f,
                                                          2.0f, 2.0f))),
                                              _mm256_set_ps(-1.0f, -1.0f, -1.0f,
                                                            -1.0f, -1.0f, -1.0f,
                                                            -1.0f, -1.0f))),
                                      _mm256_set_ps(0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(xi_6, xi_9),
                                  _mm256_set_ps(-0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f))),
                          _mm256_set_ps(omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f)),
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_add_ps(
                                                      _mm256_add_ps(
                                                          _mm256_add_ps(
                                                              _mm256_add_ps(
                                                                  _mm256_mul_ps(
                                                                      (_mm256_mul_ps(
                                                                          _mm256_add_ps(
                                                                              xi_10,
                                                                              xi_7),
                                                                          _mm256_add_ps(
                                                                              xi_10,
                                                                              xi_7))),
                                                                      _mm256_set_ps(
                                                                          0.125f,
                                                                          0.125f,
                                                                          0.125f,
                                                                          0.125f,
                                                                          0.125f,
                                                                          0.125f,
                                                                          0.125f,
                                                                          0.125f)),
                                                                  _mm256_mul_ps(
                                                                      (_mm256_mul_ps(
                                                                          xi_9,
                                                                          xi_9)),
                                                                      _mm256_set_ps(
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f))),
                                                              _mm256_mul_ps(
                                                                  xi_10,
                                                                  _mm256_set_ps(
                                                                      -0.0833333333333333f,
                                                                      -0.0833333333333333f,
                                                                      -0.0833333333333333f,
                                                                      -0.0833333333333333f,
                                                                      -0.0833333333333333f,
                                                                      -0.0833333333333333f,
                                                                      -0.0833333333333333f,
                                                                      -0.0833333333333333f))),
                                                          _mm256_mul_ps(
                                                              xi_7,
                                                              _mm256_set_ps(
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f))),
                                                      _mm256_set_ps(
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f))),
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_set_ps(
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f))),
                                          _mm256_mul_ps(
                                              xi_11,
                                              _mm256_set_ps(-1.0, -1.0, -1.0,
                                                            -1.0, -1.0, -1.0,
                                                            -1.0, -1.0))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              rho,
                                              (_mm256_mul_ps(xi_10, xi_10))),
                                          _mm256_set_ps(-0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          rho, (_mm256_mul_ps(xi_7, xi_7))),
                                      _mm256_set_ps(-0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(rho,
                                                (_mm256_mul_ps(xi_9, xi_9))),
                                  _mm256_set_ps(-0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f))),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  xi_11));
          _mm256_store_ps(
              &_data_pdfs_20_318_10[ctr_0],
              _mm256_add_ps(
                  _mm256_add_ps(
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_21,
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_mul_ps(
                                                      xi_10,
                                                      _mm256_set_ps(
                                                          2.0f, 2.0f, 2.0f,
                                                          2.0f, 2.0f, 2.0f,
                                                          2.0f, 2.0f)),
                                                  _mm256_mul_ps(
                                                      xi_7,
                                                      _mm256_set_ps(
                                                          -3.0f, -3.0f, -3.0f,
                                                          -3.0f, -3.0f, -3.0f,
                                                          -3.0f, -3.0f))),
                                              _mm256_set_ps(1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f))),
                                      _mm256_set_ps(0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f,
                                                    0.0833333333333333f)),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          xi_15,
                                          _mm256_add_ps(
                                              _mm256_add_ps(
                                                  _mm256_mul_ps(
                                                      xi_10,
                                                      _mm256_set_ps(
                                                          3.0f, 3.0f, 3.0f,
                                                          3.0f, 3.0f, 3.0f,
                                                          3.0f, 3.0f)),
                                                  _mm256_mul_ps(
                                                      xi_7,
                                                      _mm256_set_ps(
                                                          -2.0f, -2.0f, -2.0f,
                                                          -2.0f, -2.0f, -2.0f,
                                                          -2.0f, -2.0f))),
                                              _mm256_set_ps(1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f, 1.0f,
                                                            1.0f, 1.0f))),
                                      _mm256_set_ps(-0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f,
                                                    -0.0833333333333333f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(xi_6, xi_9),
                                  _mm256_set_ps(-0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f,
                                                -0.0833333333333333f))),
                          _mm256_set_ps(omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f,
                                        omega_shear * -0.5f + 1.0f)),
                      _mm256_mul_ps(
                          _mm256_add_ps(
                              _mm256_add_ps(
                                  _mm256_add_ps(
                                      _mm256_add_ps(
                                          _mm256_add_ps(
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_add_ps(
                                                      _mm256_add_ps(
                                                          _mm256_add_ps(
                                                              _mm256_add_ps(
                                                                  _mm256_mul_ps(
                                                                      (_mm256_mul_ps(
                                                                          _mm256_add_ps(
                                                                              _mm256_mul_ps(
                                                                                  xi_7,
                                                                                  _mm256_set_ps(
                                                                                      -1.0,
                                                                                      -1.0,
                                                                                      -1.0,
                                                                                      -1.0,
                                                                                      -1.0,
                                                                                      -1.0,
                                                                                      -1.0,
                                                                                      -1.0)),
                                                                              xi_10),
                                                                          _mm256_add_ps(
                                                                              _mm256_mul_ps(
                                                                                  xi_7,
                                                                                  _mm256_set_ps(
                                                                                      -1.0,
                                                                                      -1.0,
                                                                                      -1.0,
                                                                                      -1.0,
                                                                                      -1.0,
                                                                                      -1.0,
                                                                                      -1.0,
                                                                                      -1.0)),
                                                                              xi_10))),
                                                                      _mm256_set_ps(
                                                                          0.125f,
                                                                          0.125f,
                                                                          0.125f,
                                                                          0.125f,
                                                                          0.125f,
                                                                          0.125f,
                                                                          0.125f,
                                                                          0.125f)),
                                                                  _mm256_mul_ps(
                                                                      (_mm256_mul_ps(
                                                                          xi_9,
                                                                          xi_9)),
                                                                      _mm256_set_ps(
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f,
                                                                          0.0416666666666667f))),
                                                              _mm256_mul_ps(
                                                                  xi_10,
                                                                  _mm256_set_ps(
                                                                      0.0833333333333333f,
                                                                      0.0833333333333333f,
                                                                      0.0833333333333333f,
                                                                      0.0833333333333333f,
                                                                      0.0833333333333333f,
                                                                      0.0833333333333333f,
                                                                      0.0833333333333333f,
                                                                      0.0833333333333333f))),
                                                          _mm256_mul_ps(
                                                              xi_7,
                                                              _mm256_set_ps(
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f,
                                                                  -0.0833333333333333f))),
                                                      _mm256_set_ps(
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f,
                                                          -0.0138888888888889f))),
                                              _mm256_mul_ps(
                                                  rho,
                                                  _mm256_set_ps(
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f,
                                                      0.0416666666666667f))),
                                          _mm256_mul_ps(
                                              xi_2,
                                              _mm256_set_ps(-1.0, -1.0, -1.0,
                                                            -1.0, -1.0, -1.0,
                                                            -1.0, -1.0))),
                                      _mm256_mul_ps(
                                          _mm256_mul_ps(
                                              rho,
                                              (_mm256_mul_ps(xi_10, xi_10))),
                                          _mm256_set_ps(-0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f,
                                                        -0.0416666666666667f))),
                                  _mm256_mul_ps(
                                      _mm256_mul_ps(
                                          rho, (_mm256_mul_ps(xi_7, xi_7))),
                                      _mm256_set_ps(-0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f,
                                                    -0.0416666666666667f))),
                              _mm256_mul_ps(
                                  _mm256_mul_ps(rho,
                                                (_mm256_mul_ps(xi_9, xi_9))),
                                  _mm256_set_ps(-0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f,
                                                -0.0416666666666667f))),
                          _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear, omega_shear,
                                        omega_shear, omega_shear))),
                  xi_2));
        }
        for (int64_t ctr_0 = (int64_t)((_size_force_0) / (8)) * (8);
             ctr_0 < _size_force_0; ctr_0 += 1) {
          const float xi_26 = _data_pdfs_20_34_10[ctr_0];
          const float xi_27 = _data_pdfs_20_318_10[ctr_0];
          const float xi_28 = _data_pdfs_20_35_10[ctr_0];
          const float xi_29 = _data_pdfs_20_313_10[ctr_0];
          const float xi_30 = _data_pdfs_20_310_10[ctr_0];
          const float xi_31 = _data_force_20_31_10[ctr_0];
          const float xi_32 = _data_velocity_20_32_10[ctr_0];
          const float xi_33 = _data_pdfs_20_316_10[ctr_0];
          const float xi_34 = _data_velocity_20_31_10[ctr_0];
          const float xi_35 = _data_velocity_20_30_10[ctr_0];
          const float xi_36 = _data_pdfs_20_317_10[ctr_0];
          const float xi_37 = _data_pdfs_20_38_10[ctr_0];
          const float xi_38 = _data_pdfs_20_36_10[ctr_0];
          const float xi_39 = _data_pdfs_20_37_10[ctr_0];
          const float xi_40 = _data_force_20_32_10[ctr_0];
          const float xi_41 = _data_pdfs_20_312_10[ctr_0];
          const float xi_42 = _data_pdfs_20_33_10[ctr_0];
          const float xi_43 = _data_pdfs_20_30_10[ctr_0];
          const float xi_44 = _data_pdfs_20_311_10[ctr_0];
          const float xi_45 = _data_pdfs_20_315_10[ctr_0];
          const float xi_46 = _data_pdfs_20_39_10[ctr_0];
          const float xi_47 = _data_pdfs_20_314_10[ctr_0];
          const float xi_48 = _data_force_20_30_10[ctr_0];
          const float xi_49 = _data_pdfs_20_32_10[ctr_0];
          const float xi_50 = _data_pdfs_20_31_10[ctr_0];
          const float xi_1 = xi_26;
          const float xi_2 = xi_27;
          const float xi_3 = xi_28;
          const float xi_4 = xi_29;
          const float xi_5 = xi_30;
          const float xi_6 = xi_31;
          const float xi_7 = xi_32;
          const float xi_8 = xi_33;
          const float xi_9 = xi_34;
          const float xi_10 = xi_35;
          const float xi_11 = xi_36;
          const float xi_12 = xi_37;
          const float xi_13 = xi_38;
          const float xi_14 = xi_39;
          const float xi_15 = xi_40;
          const float xi_16 = xi_41;
          const float xi_17 = xi_42;
          const float xi_18 = xi_43;
          const float xi_19 = xi_45;
          const float xi_20 = xi_46;
          const float xi_21 = xi_48;
          const float xi_22 = xi_47;
          const float xi_23 = xi_44;
          const float xi_24 = xi_49;
          const float xi_25 = xi_50;
          const float rho = xi_1 + xi_11 + xi_12 + xi_13 + xi_14 + xi_16 +
                            xi_17 + xi_18 + xi_19 + xi_2 + xi_20 + xi_22 +
                            xi_23 + xi_24 + xi_25 + xi_3 + xi_4 + xi_5 + xi_8;
          const float u1Pu2 = xi_7 + xi_9;
          const float u1Mu2 = -xi_7 + xi_9;
          _data_pdfs_20_30_10[ctr_0] =
              omega_shear * (rho * (xi_10 * xi_10) * -0.333333333333333f +
                             rho * (xi_7 * xi_7) * -0.333333333333333f +
                             rho * (xi_9 * xi_9) * -0.333333333333333f +
                             rho * 0.333333333333333f - xi_18) +
              xi_18 +
              (omega_shear * -0.5f + 1.0f) *
                  (-xi_10 * xi_21 - xi_15 * xi_7 - xi_6 * xi_9);
          _data_pdfs_20_31_10[ctr_0] =
              omega_shear *
                  (rho * (xi_7 * xi_7) * -0.166666666666667f +
                   rho * (xi_9 * xi_9) * -0.166666666666667f +
                   rho *
                       ((xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                     0.050000000000000003f) *
                        (xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                     0.050000000000000003f)) *
                       -0.166666666666667f +
                   rho * ((xi_9 * xi_9) * 0.333333333333333f +
                          xi_9 * 0.166666666666667f - 0.111111111111111f) +
                   rho * 0.166666666666667f - xi_25) +
              xi_25 +
              (omega_shear * -0.5f + 1.0f) *
                  (xi_15 * xi_7 * -0.166666666666667f +
                   xi_21 *
                       (xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                    0.050000000000000003f) *
                       -0.166666666666667f +
                   xi_6 * (xi_9 * 2.0f + 1.0f) * 0.166666666666667f);
          _data_pdfs_20_32_10[ctr_0] =
              omega_shear *
                  (rho * (xi_7 * xi_7) * -0.166666666666667f +
                   rho * (xi_9 * xi_9) * -0.166666666666667f +
                   rho *
                       ((xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                     0.050000000000000003f) *
                        (xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                     0.050000000000000003f)) *
                       -0.166666666666667f +
                   rho * ((xi_9 * xi_9) * 0.333333333333333f +
                          xi_9 * -0.166666666666667f - 0.111111111111111f) +
                   rho * 0.166666666666667f - xi_24) +
              xi_24 +
              (omega_shear * -0.5f + 1.0f) *
                  (xi_15 * xi_7 * -0.166666666666667f +
                   xi_21 *
                       (xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                    0.050000000000000003f) *
                       -0.166666666666667f +
                   xi_6 * (xi_9 * 2.0f - 1.0f) * 0.166666666666667f);
          _data_pdfs_20_33_10[ctr_0] =
              omega_shear *
                  (rho * (xi_10 * xi_10) * -0.166666666666667f +
                   rho * (xi_7 * xi_7) * -0.166666666666667f +
                   rho * (xi_9 * xi_9) * -0.166666666666667f +
                   rho * ((xi_10 * xi_10) * 0.333333333333333f +
                          xi_10 * -0.166666666666667f - 0.111111111111111f) +
                   rho * 0.166666666666667f - xi_17) +
              xi_17 +
              (omega_shear * -0.5f + 1.0f) *
                  (xi_15 * xi_7 * -0.166666666666667f +
                   xi_21 * (xi_10 * 2.0f - 1.0f) * 0.166666666666667f +
                   xi_6 * xi_9 * -0.166666666666667f);
          _data_pdfs_20_34_10[ctr_0] =
              omega_shear *
                  (rho * (xi_10 * xi_10) * -0.166666666666667f +
                   rho * (xi_7 * xi_7) * -0.166666666666667f +
                   rho * (xi_9 * xi_9) * -0.166666666666667f +
                   rho * ((xi_10 * xi_10) * 0.333333333333333f +
                          xi_10 * 0.166666666666667f - 0.111111111111111f) +
                   rho * 0.166666666666667f - xi_1) +
              xi_1 +
              (omega_shear * -0.5f + 1.0f) *
                  (xi_15 * xi_7 * -0.166666666666667f +
                   xi_21 * (xi_10 * 2.0f + 1.0f) * 0.166666666666667f +
                   xi_6 * xi_9 * -0.166666666666667f);
          _data_pdfs_20_35_10[ctr_0] =
              omega_shear *
                  (rho * (xi_10 * xi_10) * -0.166666666666667f +
                   rho * (xi_7 * xi_7) * -0.166666666666667f +
                   rho * (xi_9 * xi_9) * -0.166666666666667f +
                   rho * ((xi_7 * xi_7) * 0.333333333333333f +
                          xi_7 * 0.166666666666667f - 0.111111111111111f) +
                   rho * 0.166666666666667f - xi_3) +
              xi_3 +
              (omega_shear * -0.5f + 1.0f) *
                  (xi_10 * xi_21 * -0.166666666666667f +
                   xi_15 * (xi_7 * 2.0f + 1.0f) * 0.166666666666667f +
                   xi_6 * xi_9 * -0.166666666666667f);
          _data_pdfs_20_36_10[ctr_0] =
              omega_shear *
                  (rho * (xi_10 * xi_10) * -0.166666666666667f +
                   rho * (xi_7 * xi_7) * -0.166666666666667f +
                   rho * (xi_9 * xi_9) * -0.166666666666667f +
                   rho * ((xi_7 * xi_7) * 0.333333333333333f +
                          xi_7 * -0.166666666666667f - 0.111111111111111f) +
                   rho * 0.166666666666667f - xi_13) +
              xi_13 +
              (omega_shear * -0.5f + 1.0f) *
                  (xi_10 * xi_21 * -0.166666666666667f +
                   xi_15 * (xi_7 * 2.0f - 1.0f) * 0.166666666666667f +
                   xi_6 * xi_9 * -0.166666666666667f);
          _data_pdfs_20_37_10[ctr_0] =
              omega_shear * (rho * (xi_7 * xi_7) * -0.0416666666666667f +
                             rho * (xi_9 * xi_9) * -0.0416666666666667f +
                             rho *
                                 ((xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                               0.050000000000000003f) *
                                  (xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                               0.050000000000000003f)) *
                                 -0.0416666666666667f +
                             rho * (xi_10 * -0.0833333333333333f +
                                    (xi_7 * xi_7) * 0.0416666666666667f +
                                    xi_9 * 0.0833333333333333f +
                                    ((xi_10 - xi_9 +
                                      ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                          0.050000000000000003f) *
                                     (xi_10 - xi_9 +
                                      ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                          0.050000000000000003f)) *
                                        0.125f +
                                    ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                        -0.0041666666666666666f -
                                    0.0138888888888889f) +
                             rho * 0.0416666666666667f - xi_14) +
              xi_14 +
              (omega_shear * -0.5f + 1.0f) *
                  (xi_15 * xi_7 * -0.0833333333333333f +
                   xi_21 *
                       (xi_10 * -2.0f + xi_9 * 3.0f +
                        ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                            -0.10000000000000001f +
                        1.0f) *
                       -0.0833333333333333f +
                   xi_6 *
                       (xi_10 * -3.0f + xi_9 * 2.0f +
                        ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                            -0.15000000000000002f +
                        1.0f) *
                       0.0833333333333333f);
          _data_pdfs_20_38_10[ctr_0] =
              omega_shear * (rho * (xi_7 * xi_7) * -0.0416666666666667f +
                             rho * (xi_9 * xi_9) * -0.0416666666666667f +
                             rho *
                                 ((xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                               0.050000000000000003f) *
                                  (xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                               0.050000000000000003f)) *
                                 -0.0416666666666667f +
                             rho * (xi_10 * 0.0833333333333333f +
                                    (xi_7 * xi_7) * 0.0416666666666667f +
                                    xi_9 * 0.0833333333333333f +
                                    ((xi_10 + xi_9 +
                                      ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                          0.050000000000000003f) *
                                     (xi_10 + xi_9 +
                                      ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                          0.050000000000000003f)) *
                                        0.125f +
                                    ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                        0.0041666666666666666f -
                                    0.0138888888888889f) +
                             rho * 0.0416666666666667f - xi_12) +
              xi_12 +
              (omega_shear * -0.5f + 1.0f) *
                  (xi_15 * xi_7 * -0.0833333333333333f +
                   xi_21 *
                       (xi_10 * 2.0f + xi_9 * 3.0f +
                        ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                            0.10000000000000001f +
                        1.0f) *
                       0.0833333333333333f +
                   xi_6 *
                       (xi_10 * 3.0f + xi_9 * 2.0f +
                        ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                            0.15000000000000002f +
                        1.0f) *
                       0.0833333333333333f);
          _data_pdfs_20_39_10[ctr_0] =
              omega_shear * (rho * (xi_7 * xi_7) * -0.0416666666666667f +
                             rho * (xi_9 * xi_9) * -0.0416666666666667f +
                             rho *
                                 ((xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                               0.050000000000000003f) *
                                  (xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                               0.050000000000000003f)) *
                                 -0.0416666666666667f +
                             rho * (xi_10 * -0.0833333333333333f +
                                    (xi_7 * xi_7) * 0.0416666666666667f +
                                    xi_9 * -0.0833333333333333f +
                                    ((xi_10 + xi_9 +
                                      ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                          0.050000000000000003f) *
                                     (xi_10 + xi_9 +
                                      ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                          0.050000000000000003f)) *
                                        0.125f +
                                    ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                        -0.0041666666666666666f -
                                    0.0138888888888889f) +
                             rho * 0.0416666666666667f - xi_20) +
              xi_20 +
              (omega_shear * -0.5f + 1.0f) *
                  (xi_15 * xi_7 * -0.0833333333333333f +
                   xi_21 *
                       (xi_10 * 2.0f + xi_9 * 3.0f +
                        ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                            0.10000000000000001f -
                        1.0f) *
                       0.0833333333333333f +
                   xi_6 *
                       (xi_10 * 3.0f + xi_9 * 2.0f +
                        ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                            0.15000000000000002f -
                        1.0f) *
                       0.0833333333333333f);
          _data_pdfs_20_310_10[ctr_0] =
              omega_shear * (rho * (xi_7 * xi_7) * -0.0416666666666667f +
                             rho * (xi_9 * xi_9) * -0.0416666666666667f +
                             rho *
                                 ((xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                               0.050000000000000003f) *
                                  (xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                               0.050000000000000003f)) *
                                 -0.0416666666666667f +
                             rho * (xi_10 * 0.0833333333333333f +
                                    (xi_7 * xi_7) * 0.0416666666666667f +
                                    xi_9 * -0.0833333333333333f +
                                    ((xi_10 - xi_9 +
                                      ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                          0.050000000000000003f) *
                                     (xi_10 - xi_9 +
                                      ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                          0.050000000000000003f)) *
                                        0.125f +
                                    ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                        0.0041666666666666666f -
                                    0.0138888888888889f) +
                             rho * 0.0416666666666667f - xi_5) +
              xi_5 +
              (omega_shear * -0.5f + 1.0f) *
                  (xi_15 * xi_7 * -0.0833333333333333f +
                   xi_21 *
                       (xi_10 * 2.0f + xi_9 * -3.0f +
                        ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                            0.10000000000000001f +
                        1.0f) *
                       0.0833333333333333f +
                   xi_6 *
                       (xi_10 * 3.0f + xi_9 * -2.0f +
                        ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                            0.15000000000000002f +
                        1.0f) *
                       -0.0833333333333333f);
          _data_pdfs_20_311_10[ctr_0] =
              omega_shear *
                  (rho * (xi_7 * xi_7) * -0.0416666666666667f +
                   rho * (xi_9 * xi_9) * -0.0416666666666667f +
                   rho *
                       ((xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                     0.050000000000000003f) *
                        (xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                     0.050000000000000003f)) *
                       -0.0416666666666667f +
                   rho *
                       ((u1Pu2 * u1Pu2) * 0.125f + u1Pu2 * 0.0833333333333333f +
                        ((xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                      0.050000000000000003f) *
                         (xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                      0.050000000000000003f)) *
                            0.0416666666666667f -
                        0.0138888888888889f) +
                   rho * 0.0416666666666667f - xi_23) +
              xi_23 +
              (omega_shear * -0.5f + 1.0f) *
                  (xi_15 * (xi_7 * 2.0f + xi_9 * 3.0f + 1.0f) *
                       0.0833333333333333f +
                   xi_21 *
                       (xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                    0.050000000000000003f) *
                       -0.0833333333333333f +
                   xi_6 * (xi_7 * 3.0f + xi_9 * 2.0f + 1.0f) *
                       0.0833333333333333f);
          _data_pdfs_20_312_10[ctr_0] =
              omega_shear * (rho * (xi_7 * xi_7) * -0.0416666666666667f +
                             rho * (xi_9 * xi_9) * -0.0416666666666667f +
                             rho *
                                 ((xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                               0.050000000000000003f) *
                                  (xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                               0.050000000000000003f)) *
                                 -0.0416666666666667f +
                             rho * ((u1Mu2 * u1Mu2) * 0.125f +
                                    u1Mu2 * -0.0833333333333333f +
                                    ((xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                                  0.050000000000000003f) *
                                     (xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                                  0.050000000000000003f)) *
                                        0.0416666666666667f -
                                    0.0138888888888889f) +
                             rho * 0.0416666666666667f - xi_16) +
              xi_16 +
              (omega_shear * -0.5f + 1.0f) *
                  (xi_15 * (xi_7 * 2.0f + xi_9 * -3.0f + 1.0f) *
                       0.0833333333333333f +
                   xi_21 *
                       (xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                    0.050000000000000003f) *
                       -0.0833333333333333f +
                   xi_6 * (xi_7 * 3.0f + xi_9 * -2.0f + 1.0f) *
                       -0.0833333333333333f);
          _data_pdfs_20_313_10[ctr_0] =
              omega_shear * (rho * (xi_10 * xi_10) * -0.0416666666666667f +
                             rho * (xi_7 * xi_7) * -0.0416666666666667f +
                             rho * (xi_9 * xi_9) * -0.0416666666666667f +
                             rho * (xi_10 * -0.0833333333333333f +
                                    xi_7 * 0.0833333333333333f +
                                    (xi_9 * xi_9) * 0.0416666666666667f +
                                    ((xi_10 - xi_7) * (xi_10 - xi_7)) * 0.125f -
                                    0.0138888888888889f) +
                             rho * 0.0416666666666667f - xi_4) +
              xi_4 +
              (omega_shear * -0.5f + 1.0f) *
                  (xi_15 * (xi_10 * -3.0f + xi_7 * 2.0f + 1.0f) *
                       0.0833333333333333f +
                   xi_21 * (xi_10 * -2.0f + xi_7 * 3.0f + 1.0f) *
                       -0.0833333333333333f +
                   xi_6 * xi_9 * -0.0833333333333333f);
          _data_pdfs_20_314_10[ctr_0] =
              omega_shear * (rho * (xi_10 * xi_10) * -0.0416666666666667f +
                             rho * (xi_7 * xi_7) * -0.0416666666666667f +
                             rho * (xi_9 * xi_9) * -0.0416666666666667f +
                             rho * (xi_10 * 0.0833333333333333f +
                                    xi_7 * 0.0833333333333333f +
                                    (xi_9 * xi_9) * 0.0416666666666667f +
                                    ((xi_10 + xi_7) * (xi_10 + xi_7)) * 0.125f -
                                    0.0138888888888889f) +
                             rho * 0.0416666666666667f - xi_22) +
              xi_22 +
              (omega_shear * -0.5f + 1.0f) *
                  (xi_15 * (xi_10 * 3.0f + xi_7 * 2.0f + 1.0f) *
                       0.0833333333333333f +
                   xi_21 * (xi_10 * 2.0f + xi_7 * 3.0f + 1.0f) *
                       0.0833333333333333f +
                   xi_6 * xi_9 * -0.0833333333333333f);
          _data_pdfs_20_315_10[ctr_0] =
              omega_shear *
                  (rho * (xi_7 * xi_7) * -0.0416666666666667f +
                   rho * (xi_9 * xi_9) * -0.0416666666666667f +
                   rho *
                       ((xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                     0.050000000000000003f) *
                        (xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                     0.050000000000000003f)) *
                       -0.0416666666666667f +
                   rho *
                       ((u1Mu2 * u1Mu2) * 0.125f + u1Mu2 * 0.0833333333333333f +
                        ((xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                      0.050000000000000003f) *
                         (xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                      0.050000000000000003f)) *
                            0.0416666666666667f -
                        0.0138888888888889f) +
                   rho * 0.0416666666666667f - xi_19) +
              xi_19 +
              (omega_shear * -0.5f + 1.0f) *
                  (xi_15 * (xi_7 * -2.0f + xi_9 * 3.0f + 1.0f) *
                       -0.0833333333333333f +
                   xi_21 *
                       (xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                    0.050000000000000003f) *
                       -0.0833333333333333f +
                   xi_6 * (xi_7 * -3.0f + xi_9 * 2.0f + 1.0f) *
                       0.0833333333333333f);
          _data_pdfs_20_316_10[ctr_0] =
              omega_shear * (rho * (xi_7 * xi_7) * -0.0416666666666667f +
                             rho * (xi_9 * xi_9) * -0.0416666666666667f +
                             rho *
                                 ((xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                               0.050000000000000003f) *
                                  (xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                               0.050000000000000003f)) *
                                 -0.0416666666666667f +
                             rho * ((u1Pu2 * u1Pu2) * 0.125f +
                                    u1Pu2 * -0.0833333333333333f +
                                    ((xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                                  0.050000000000000003f) *
                                     (xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                                  0.050000000000000003f)) *
                                        0.0416666666666667f -
                                    0.0138888888888889f) +
                             rho * 0.0416666666666667f - xi_8) +
              xi_8 +
              (omega_shear * -0.5f + 1.0f) *
                  (xi_15 * (xi_7 * 2.0f + xi_9 * 3.0f - 1.0f) *
                       0.0833333333333333f +
                   xi_21 *
                       (xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                    0.050000000000000003f) *
                       -0.0833333333333333f +
                   xi_6 * (xi_7 * 3.0f + xi_9 * 2.0f - 1.0f) *
                       0.0833333333333333f);
          _data_pdfs_20_317_10[ctr_0] =
              omega_shear * (rho * (xi_10 * xi_10) * -0.0416666666666667f +
                             rho * (xi_7 * xi_7) * -0.0416666666666667f +
                             rho * (xi_9 * xi_9) * -0.0416666666666667f +
                             rho * (xi_10 * -0.0833333333333333f +
                                    xi_7 * -0.0833333333333333f +
                                    (xi_9 * xi_9) * 0.0416666666666667f +
                                    ((xi_10 + xi_7) * (xi_10 + xi_7)) * 0.125f -
                                    0.0138888888888889f) +
                             rho * 0.0416666666666667f - xi_11) +
              xi_11 +
              (omega_shear * -0.5f + 1.0f) *
                  (xi_15 * (xi_10 * 3.0f + xi_7 * 2.0f - 1.0f) *
                       0.0833333333333333f +
                   xi_21 * (xi_10 * 2.0f + xi_7 * 3.0f - 1.0f) *
                       0.0833333333333333f +
                   xi_6 * xi_9 * -0.0833333333333333f);
          _data_pdfs_20_318_10[ctr_0] =
              omega_shear * (rho * (xi_10 * xi_10) * -0.0416666666666667f +
                             rho * (xi_7 * xi_7) * -0.0416666666666667f +
                             rho * (xi_9 * xi_9) * -0.0416666666666667f +
                             rho * (xi_10 * 0.0833333333333333f +
                                    xi_7 * -0.0833333333333333f +
                                    (xi_9 * xi_9) * 0.0416666666666667f +
                                    ((xi_10 - xi_7) * (xi_10 - xi_7)) * 0.125f -
                                    0.0138888888888889f) +
                             rho * 0.0416666666666667f - xi_2) +
              xi_2 +
              (omega_shear * -0.5f + 1.0f) *
                  (xi_15 * (xi_10 * 3.0f + xi_7 * -2.0f + 1.0f) *
                       -0.0833333333333333f +
                   xi_21 * (xi_10 * 2.0f + xi_7 * -3.0f + 1.0f) *
                       0.0833333333333333f +
                   xi_6 * xi_9 * -0.0833333333333333f);
        }
      }
    }
  }
}
} // namespace internal_9a18f2f4073cdcc5365cdfddb752069e

void CollideSweepSinglePrecisionLeesEdwardsAVX::run(IBlock *block) {
  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);
  auto force = block->getData<field::GhostLayerField<float, 3>>(forceID);
  auto velocity = block->getData<field::GhostLayerField<float, 3>>(velocityID);

  auto &omega_shear = this->omega_shear_;
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force->nrOfGhostLayers()));
  float *RESTRICT const _data_force = force->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()));
  float *RESTRICT _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(velocity->nrOfGhostLayers()));
  float *RESTRICT const _data_velocity = velocity->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx);
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
  const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
  const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
  const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
  internal_9a18f2f4073cdcc5365cdfddb752069e::
      collidesweepsingleprecisionleesedwardsavx_collidesweepsingleprecisionleesedwardsavx(
          _data_force, _data_pdfs, _data_velocity, _size_force_0, _size_force_1,
          _size_force_2, _stride_force_1, _stride_force_2, _stride_force_3,
          _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_velocity_1,
          _stride_velocity_2, _stride_velocity_3, omega_shear);
}

void CollideSweepSinglePrecisionLeesEdwardsAVX::runOnCellInterval(
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

  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);
  auto force = block->getData<field::GhostLayerField<float, 3>>(forceID);
  auto velocity = block->getData<field::GhostLayerField<float, 3>>(velocityID);

  auto &omega_shear = this->omega_shear_;
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(force->nrOfGhostLayers()));
  float *RESTRICT const _data_force =
      force->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
  float *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(velocity->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(velocity->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(velocity->nrOfGhostLayers()));
  float *RESTRICT const _data_velocity =
      velocity->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx);
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
  const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
  const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
  const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
  internal_9a18f2f4073cdcc5365cdfddb752069e::
      collidesweepsingleprecisionleesedwardsavx_collidesweepsingleprecisionleesedwardsavx(
          _data_force, _data_pdfs, _data_velocity, _size_force_0, _size_force_1,
          _size_force_2, _stride_force_1, _stride_force_2, _stride_force_3,
          _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_velocity_1,
          _stride_velocity_2, _stride_velocity_3, omega_shear);
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