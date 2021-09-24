// kernel generated with pystencils v0.3.3+44.gf2d4181, lbmpy
// v0.3.1+53.g2faceda, lbmpy_walberla/pystencils_walberla from commit
// b17ca5caf00db7d19f86c5f85c6f67fec6c16aff

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
//! \\file CollideSweepLeesEdwardsAVX.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

#include <cmath>

#include "CollideSweepLeesEdwardsAVX.h"
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

namespace internal_collidesweepleesedwardsavx {
static FUNC_PREFIX void collidesweepleesedwardsavx(
    double *RESTRICT const _data_density, double *RESTRICT const _data_force,
    double *RESTRICT _data_pdfs, double *RESTRICT const _data_velocity,
    int64_t const _size_density_0, int64_t const _size_density_1,
    int64_t const _size_density_2, int64_t const _stride_density_1,
    int64_t const _stride_density_2, int64_t const _stride_force_1,
    int64_t const _stride_force_2, int64_t const _stride_force_3,
    int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
    int64_t const _stride_pdfs_3, int64_t const _stride_velocity_1,
    int64_t const _stride_velocity_2, int64_t const _stride_velocity_3,
    double omega_bulk, double omega_even, double omega_odd, double omega_shear,
    bool points_down, bool points_up) {
  const double xi_10 = -omega_shear + 2.0;
  const double xi_11 = xi_10 * 0.5;
  const double xi_16 = xi_10 * 0.0833333333333333;
  const double xi_21 = xi_10 * 0.166666666666667;
  const double xi_31 = xi_10 * 0.25;
  const double xi_36 = xi_10 * 0.0416666666666667;
  const double xi_103 = omega_odd * 0.25;
  const double xi_110 = omega_odd * 0.0833333333333333;
  const double xi_151 = omega_shear * 0.25;
  const double xi_173 = omega_odd * 0.0416666666666667;
  const double xi_175 = omega_odd * 0.125;
  const int64_t rr_0 = 0.0;
  const double xi_116 = rr_0 * 0.166666666666667;
  const double xi_156 = rr_0 * 0.0833333333333333;
  for (int64_t ctr_2 = 0; ctr_2 < _size_density_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_20_312 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_313 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_36 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_35 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    double *RESTRICT _data_density_20_30 =
        _data_density + _stride_density_2 * ctr_2;
    double *RESTRICT _data_velocity_20_32 =
        _data_velocity + _stride_velocity_2 * ctr_2 + 2 * _stride_velocity_3;
    double *RESTRICT _data_pdfs_20_316 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    double *RESTRICT _data_velocity_20_30 =
        _data_velocity + _stride_velocity_2 * ctr_2;
    double *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    double *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_318 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_315 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_311 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3;
    double *RESTRICT _data_force_20_31 =
        _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    double *RESTRICT _data_pdfs_20_317 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_33 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_32 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    double *RESTRICT _data_force_20_32 =
        _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    double *RESTRICT _data_pdfs_20_31 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_34 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_314 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3;
    double *RESTRICT _data_velocity_20_31 =
        _data_velocity + _stride_velocity_2 * ctr_2 + _stride_velocity_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_density_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_20_312_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_312;
      double *RESTRICT _data_pdfs_20_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_313;
      double *RESTRICT _data_pdfs_20_36_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_36;
      double *RESTRICT _data_pdfs_20_35_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_35;
      double *RESTRICT _data_pdfs_20_310_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_310;
      double *RESTRICT _data_pdfs_20_30_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      double *RESTRICT _data_density_20_30_10 =
          _stride_density_1 * ctr_1 + _data_density_20_30;
      double *RESTRICT _data_velocity_20_32_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_32;
      double *RESTRICT _data_pdfs_20_316_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_316;
      double *RESTRICT _data_pdfs_20_38_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_38;
      double *RESTRICT _data_velocity_20_30_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_30;
      double *RESTRICT _data_force_20_30_10 =
          _stride_force_1 * ctr_1 + _data_force_20_30;
      double *RESTRICT _data_pdfs_20_37_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_37;
      double *RESTRICT _data_pdfs_20_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_318;
      double *RESTRICT _data_pdfs_20_315_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_315;
      double *RESTRICT _data_pdfs_20_311_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_311;
      double *RESTRICT _data_force_20_31_10 =
          _stride_force_1 * ctr_1 + _data_force_20_31;
      double *RESTRICT _data_pdfs_20_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_317;
      double *RESTRICT _data_pdfs_20_33_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      double *RESTRICT _data_pdfs_20_32_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_32;
      double *RESTRICT _data_force_20_32_10 =
          _stride_force_1 * ctr_1 + _data_force_20_32;
      double *RESTRICT _data_pdfs_20_31_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_31;
      double *RESTRICT _data_pdfs_20_34_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      double *RESTRICT _data_pdfs_20_39_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_39;
      double *RESTRICT _data_pdfs_20_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_314;
      double *RESTRICT _data_velocity_20_31_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_31;
      for (int64_t ctr_0 = 0;
           ctr_0 < ((_size_density_0) % (4) == 0
                        ? _size_density_0
                        : ((int64_t)((_size_density_0) / (4)) + 1) * (4));
           ctr_0 += 4) {
        const __m256d xi_199 = _mm256_load_pd(&_data_pdfs_20_312_10[ctr_0]);
        const __m256d xi_200 = _mm256_load_pd(&_data_pdfs_20_313_10[ctr_0]);
        const __m256d xi_201 = _mm256_load_pd(&_data_pdfs_20_36_10[ctr_0]);
        const __m256d xi_202 = _mm256_load_pd(&_data_pdfs_20_35_10[ctr_0]);
        const __m256d xi_203 = _mm256_load_pd(&_data_pdfs_20_310_10[ctr_0]);
        const __m256d xi_204 = _mm256_load_pd(&_data_pdfs_20_30_10[ctr_0]);
        const __m256d xi_205 = _mm256_load_pd(&_data_density_20_30_10[ctr_0]);
        const __m256d xi_206 = _mm256_load_pd(&_data_velocity_20_32_10[ctr_0]);
        const __m256d xi_207 = _mm256_load_pd(&_data_pdfs_20_316_10[ctr_0]);
        const __m256d xi_208 = _mm256_load_pd(&_data_pdfs_20_38_10[ctr_0]);
        const __m256d xi_209 = _mm256_load_pd(&_data_velocity_20_30_10[ctr_0]);
        const __m256d xi_210 = _mm256_load_pd(&_data_force_20_30_10[ctr_0]);
        const __m256d xi_211 = _mm256_load_pd(&_data_pdfs_20_37_10[ctr_0]);
        const __m256d xi_212 = _mm256_load_pd(&_data_pdfs_20_318_10[ctr_0]);
        const __m256d xi_213 = _mm256_load_pd(&_data_pdfs_20_315_10[ctr_0]);
        const __m256d xi_214 = _mm256_load_pd(&_data_pdfs_20_311_10[ctr_0]);
        const __m256d xi_215 = _mm256_load_pd(&_data_force_20_31_10[ctr_0]);
        const __m256d xi_216 = _mm256_load_pd(&_data_pdfs_20_317_10[ctr_0]);
        const __m256d xi_217 = _mm256_load_pd(&_data_pdfs_20_33_10[ctr_0]);
        const __m256d xi_218 = _mm256_load_pd(&_data_pdfs_20_32_10[ctr_0]);
        const __m256d xi_219 = _mm256_load_pd(&_data_force_20_32_10[ctr_0]);
        const __m256d xi_220 = _mm256_load_pd(&_data_pdfs_20_31_10[ctr_0]);
        const __m256d xi_221 = _mm256_load_pd(&_data_pdfs_20_34_10[ctr_0]);
        const __m256d xi_222 = _mm256_load_pd(&_data_pdfs_20_39_10[ctr_0]);
        const __m256d xi_223 = _mm256_load_pd(&_data_pdfs_20_314_10[ctr_0]);
        const __m256d xi_224 = _mm256_load_pd(&_data_velocity_20_31_10[ctr_0]);
        const __m256d xi_15 = _mm256_mul_pd(
            xi_215, _mm256_set_pd(0.166666666666667, 0.166666666666667,
                                  0.166666666666667, 0.166666666666667));
        const __m256d xi_23 = _mm256_mul_pd(
            xi_210, _mm256_set_pd(0.166666666666667, 0.166666666666667,
                                  0.166666666666667, 0.166666666666667));
        const __m256d xi_27 = _mm256_mul_pd(
            xi_219, _mm256_set_pd(0.166666666666667, 0.166666666666667,
                                  0.166666666666667, 0.166666666666667));
        const __m256d xi_30 =
            _mm256_mul_pd(xi_215, _mm256_set_pd(0.5, 0.5, 0.5, 0.5));
        const __m256d xi_34 = _mm256_mul_pd(
            xi_210, _mm256_set_pd(0.0833333333333333, 0.0833333333333333,
                                  0.0833333333333333, 0.0833333333333333));
        const __m256d xi_38 = _mm256_mul_pd(
            xi_215, _mm256_set_pd(0.0833333333333333, 0.0833333333333333,
                                  0.0833333333333333, 0.0833333333333333));
        const __m256d xi_48 = _mm256_mul_pd(
            xi_219, _mm256_set_pd(0.0833333333333333, 0.0833333333333333,
                                  0.0833333333333333, 0.0833333333333333));
        const __m256d xi_59 =
            _mm256_mul_pd(xi_204, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_60 = _mm256_add_pd(
            _mm256_add_pd(
                _mm256_mul_pd(xi_201, _mm256_set_pd(3.0, 3.0, 3.0, 3.0)),
                _mm256_mul_pd(xi_202, _mm256_set_pd(3.0, 3.0, 3.0, 3.0))),
            xi_59);
        const __m256d xi_61 = _mm256_mul_pd(
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_add_pd(
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_mul_pd(
                                        xi_199,
                                        _mm256_set_pd(-3.0, -3.0, -3.0, -3.0)),
                                    _mm256_mul_pd(
                                        xi_207,
                                        _mm256_set_pd(-3.0, -3.0, -3.0, -3.0))),
                                _mm256_mul_pd(
                                    xi_213,
                                    _mm256_set_pd(-3.0, -3.0, -3.0, -3.0))),
                            _mm256_mul_pd(
                                xi_214, _mm256_set_pd(-3.0, -3.0, -3.0, -3.0))),
                        _mm256_mul_pd(xi_218,
                                      _mm256_set_pd(3.0, 3.0, 3.0, 3.0))),
                    _mm256_mul_pd(xi_220, _mm256_set_pd(3.0, 3.0, 3.0, 3.0))),
                xi_60),
            _mm256_set_pd(omega_even, omega_even, omega_even, omega_even));
        const __m256d xi_62 = _mm256_add_pd(
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_mul_pd(xi_199, _mm256_set_pd(2.0, 2.0, 2.0, 2.0)),
                    _mm256_mul_pd(xi_207, _mm256_set_pd(2.0, 2.0, 2.0, 2.0))),
                _mm256_mul_pd(xi_213, _mm256_set_pd(2.0, 2.0, 2.0, 2.0))),
            _mm256_mul_pd(xi_214, _mm256_set_pd(2.0, 2.0, 2.0, 2.0)));
        const __m256d xi_63 = _mm256_add_pd(
            _mm256_add_pd(
                _mm256_mul_pd(xi_217, _mm256_set_pd(5.0, 5.0, 5.0, 5.0)),
                _mm256_mul_pd(xi_221, _mm256_set_pd(5.0, 5.0, 5.0, 5.0))),
            xi_62);
        const __m256d xi_64 = _mm256_mul_pd(
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_add_pd(
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_add_pd(
                                        _mm256_mul_pd(
                                            xi_200, _mm256_set_pd(-5.0, -5.0,
                                                                  -5.0, -5.0)),
                                        _mm256_mul_pd(
                                            xi_212, _mm256_set_pd(-5.0, -5.0,
                                                                  -5.0, -5.0))),
                                    _mm256_mul_pd(
                                        xi_216,
                                        _mm256_set_pd(-5.0, -5.0, -5.0, -5.0))),
                                _mm256_mul_pd(
                                    xi_218,
                                    _mm256_set_pd(-2.0, -2.0, -2.0, -2.0))),
                            _mm256_mul_pd(
                                xi_220, _mm256_set_pd(-2.0, -2.0, -2.0, -2.0))),
                        _mm256_mul_pd(xi_223,
                                      _mm256_set_pd(-5.0, -5.0, -5.0, -5.0))),
                    xi_60),
                xi_63),
            _mm256_set_pd(omega_even, omega_even, omega_even, omega_even));
        const __m256d xi_67 =
            _mm256_mul_pd(xi_222, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_69 =
            _mm256_mul_pd(xi_208, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_70 = _mm256_add_pd(
            _mm256_mul_pd(xi_203, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
            xi_69);
        const __m256d xi_72 =
            _mm256_mul_pd(xi_216, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_73 =
            _mm256_mul_pd(xi_212, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_74 =
            _mm256_mul_pd(xi_214, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_75 =
            _mm256_mul_pd(xi_199, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_76 = _mm256_add_pd(xi_74, xi_75);
        const __m256d xi_77 =
            _mm256_mul_pd(xi_207, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_78 =
            _mm256_mul_pd(xi_213, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_79 = _mm256_add_pd(xi_77, xi_78);
        const __m256d xi_80 = _mm256_add_pd(
            _mm256_add_pd(_mm256_add_pd(xi_72, xi_73), xi_76), xi_79);
        const __m256d xi_81 =
            _mm256_mul_pd(xi_211, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_82 =
            _mm256_mul_pd(xi_223, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_83 =
            _mm256_mul_pd(xi_200, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_84 = _mm256_add_pd(xi_82, xi_83);
        const __m256d xi_85 = _mm256_add_pd(xi_81, xi_84);
        const __m256d xi_87 =
            _mm256_mul_pd(xi_200, _mm256_set_pd(2.0, 2.0, 2.0, 2.0));
        const __m256d xi_88 =
            _mm256_mul_pd(xi_223, _mm256_set_pd(2.0, 2.0, 2.0, 2.0));
        const __m256d xi_89 = _mm256_add_pd(
            _mm256_mul_pd(xi_212, _mm256_set_pd(2.0, 2.0, 2.0, 2.0)),
            _mm256_mul_pd(xi_216, _mm256_set_pd(2.0, 2.0, 2.0, 2.0)));
        const __m256d xi_90 = _mm256_mul_pd(
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
                                                                xi_201,
                                                                _mm256_set_pd(
                                                                    -4.0, -4.0,
                                                                    -4.0,
                                                                    -4.0)),
                                                            _mm256_mul_pd(
                                                                xi_202,
                                                                _mm256_set_pd(
                                                                    -4.0, -4.0,
                                                                    -4.0,
                                                                    -4.0))),
                                                        _mm256_mul_pd(
                                                            xi_203,
                                                            _mm256_set_pd(
                                                                -7.0, -7.0,
                                                                -7.0, -7.0))),
                                                    _mm256_mul_pd(
                                                        xi_208,
                                                        _mm256_set_pd(
                                                            -7.0, -7.0, -7.0,
                                                            -7.0))),
                                                _mm256_mul_pd(
                                                    xi_211,
                                                    _mm256_set_pd(-7.0, -7.0,
                                                                  -7.0, -7.0))),
                                            _mm256_mul_pd(
                                                xi_218,
                                                _mm256_set_pd(5.0, 5.0, 5.0,
                                                              5.0))),
                                        _mm256_mul_pd(
                                            xi_220,
                                            _mm256_set_pd(5.0, 5.0, 5.0, 5.0))),
                                    _mm256_mul_pd(
                                        xi_222,
                                        _mm256_set_pd(-7.0, -7.0, -7.0, -7.0))),
                                xi_59),
                            xi_63),
                        xi_87),
                    xi_88),
                xi_89),
            _mm256_set_pd(omega_even, omega_even, omega_even, omega_even));
        const __m256d xi_91 =
            _mm256_mul_pd(xi_217, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_92 =
            _mm256_mul_pd(xi_201, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_93 = _mm256_add_pd(
            _mm256_mul_pd(xi_221, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
            xi_70);
        const __m256d xi_94 = _mm256_add_pd(
            _mm256_mul_pd(xi_202, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
            xi_84);
        const __m256d xi_95 =
            _mm256_mul_pd(xi_218, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_96 = _mm256_add_pd(
            _mm256_mul_pd(xi_220, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
            xi_81);
        const __m256d xi_97 = _mm256_add_pd(xi_95, xi_96);
        const __m256d xi_99 = _mm256_add_pd(xi_220, xi_95);
        const __m256d xi_100 = _mm256_add_pd(xi_199, xi_74);
        const __m256d xi_101 =
            _mm256_add_pd(_mm256_add_pd(xi_100, xi_207), xi_78);
        const __m256d xi_102 = _mm256_add_pd(xi_101, xi_99);
        const __m256d xi_104 = _mm256_mul_pd(
            xi_102, _mm256_set_pd(xi_103, xi_103, xi_103, xi_103));
        const __m256d xi_105 =
            _mm256_mul_pd(xi_211, _mm256_set_pd(2.0, 2.0, 2.0, 2.0));
        const __m256d xi_106 =
            _mm256_mul_pd(xi_203, _mm256_set_pd(2.0, 2.0, 2.0, 2.0));
        const __m256d xi_107 = _mm256_add_pd(xi_213, xi_77);
        const __m256d xi_108 = _mm256_add_pd(
            _mm256_mul_pd(xi_208, _mm256_set_pd(-2.0, -2.0, -2.0, -2.0)),
            _mm256_mul_pd(xi_222, _mm256_set_pd(2.0, 2.0, 2.0, 2.0)));
        const __m256d xi_109 = _mm256_add_pd(
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_add_pd(
                            _mm256_add_pd(
                                _mm256_mul_pd(
                                    xi_105,
                                    _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                                xi_106),
                            xi_107),
                        xi_108),
                    xi_214),
                xi_75),
            xi_99);
        const __m256d xi_111 = _mm256_mul_pd(
            xi_109, _mm256_set_pd(xi_110, xi_110, xi_110, xi_110));
        const __m256d xi_112 =
            _mm256_mul_pd(xi_111, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_114 = _mm256_add_pd(xi_203, xi_69);
        const __m256d xi_118 =
            _mm256_add_pd(_mm256_add_pd(xi_201, xi_212), xi_216);
        const __m256d xi_129 = _mm256_add_pd(xi_200, xi_82);
        const __m256d xi_130 =
            _mm256_add_pd(_mm256_add_pd(xi_129, xi_216), xi_73);
        const __m256d xi_133 = _mm256_add_pd(xi_221, xi_91);
        const __m256d xi_134 = _mm256_add_pd(xi_130, xi_133);
        const __m256d xi_135 = _mm256_mul_pd(
            xi_134, _mm256_set_pd(xi_103, xi_103, xi_103, xi_103));
        const __m256d xi_136 = _mm256_add_pd(xi_212, xi_72);
        const __m256d xi_137 = _mm256_add_pd(
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_add_pd(
                            _mm256_add_pd(
                                _mm256_mul_pd(
                                    xi_106,
                                    _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                                xi_105),
                            xi_108),
                        xi_133),
                    xi_136),
                xi_223),
            xi_83);
        const __m256d xi_138 = _mm256_mul_pd(
            xi_137, _mm256_set_pd(xi_110, xi_110, xi_110, xi_110));
        const __m256d xi_140 =
            _mm256_mul_pd(xi_138, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_141 = _mm256_add_pd(xi_202, xi_92);
        const __m256d xi_142 =
            _mm256_add_pd(_mm256_add_pd(xi_207, xi_213), xi_76);
        const __m256d xi_143 = _mm256_add_pd(xi_141, xi_142);
        const __m256d xi_144 = _mm256_mul_pd(
            xi_143, _mm256_set_pd(xi_103, xi_103, xi_103, xi_103));
        const __m256d xi_145 = _mm256_add_pd(
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_add_pd(
                            _mm256_add_pd(
                                _mm256_mul_pd(xi_87, _mm256_set_pd(-1.0, -1.0,
                                                                   -1.0, -1.0)),
                                _mm256_mul_pd(
                                    xi_88,
                                    _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))),
                            xi_141),
                        xi_199),
                    xi_214),
                xi_79),
            xi_89);
        const __m256d xi_146 = _mm256_mul_pd(
            xi_145, _mm256_set_pd(xi_110, xi_110, xi_110, xi_110));
        const __m256d xi_147 =
            _mm256_mul_pd(xi_146, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_154 = _mm256_mul_pd(
            xi_90, _mm256_set_pd(0.0138888888888889, 0.0138888888888889,
                                 0.0138888888888889, 0.0138888888888889));
        const __m256d xi_169 = _mm256_mul_pd(
            xi_64, _mm256_set_pd(-0.00714285714285714, -0.00714285714285714,
                                 -0.00714285714285714, -0.00714285714285714));
        const __m256d xi_171 =
            _mm256_mul_pd(xi_61, _mm256_set_pd(0.025, 0.025, 0.025, 0.025));
        const __m256d xi_174 = _mm256_mul_pd(
            xi_145, _mm256_set_pd(xi_173, xi_173, xi_173, xi_173));
        const __m256d xi_176 = _mm256_mul_pd(
            xi_143, _mm256_set_pd(xi_175, xi_175, xi_175, xi_175));
        const __m256d xi_177 = _mm256_mul_pd(
            xi_90, _mm256_set_pd(-0.00396825396825397, -0.00396825396825397,
                                 -0.00396825396825397, -0.00396825396825397));
        const __m256d xi_181 = _mm256_mul_pd(
            xi_109, _mm256_set_pd(xi_173, xi_173, xi_173, xi_173));
        const __m256d xi_182 = _mm256_mul_pd(
            xi_102, _mm256_set_pd(xi_175, xi_175, xi_175, xi_175));
        const __m256d xi_188 = _mm256_mul_pd(
            xi_64, _mm256_set_pd(0.0178571428571429, 0.0178571428571429,
                                 0.0178571428571429, 0.0178571428571429));
        const __m256d xi_191 = _mm256_mul_pd(
            xi_134, _mm256_set_pd(xi_175, xi_175, xi_175, xi_175));
        const __m256d xi_192 = _mm256_mul_pd(
            xi_137, _mm256_set_pd(xi_173, xi_173, xi_173, xi_173));
        const __m256d rho = xi_205;
        const __m256d xi_98 = _mm256_mul_pd(
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_add_pd(
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_add_pd(_mm256_add_pd(rho, xi_59),
                                                  xi_67),
                                    xi_80),
                                xi_91),
                            xi_92),
                        xi_93),
                    xi_94),
                xi_97),
            _mm256_set_epi32(rr_0, rr_0, rr_0, rr_0, rr_0, rr_0, rr_0, rr_0));
        const __m256d xi_122 = _mm256_mul_pd(
            xi_98, _mm256_set_pd(0.0555555555555556, 0.0555555555555556,
                                 0.0555555555555556, 0.0555555555555556));
        const __m256d xi_123 = _mm256_add_pd(
            _mm256_mul_pd(
                xi_90, _mm256_set_pd(-0.0198412698412698, -0.0198412698412698,
                                     -0.0198412698412698, -0.0198412698412698)),
            xi_122);
        const __m256d u_0 = _mm256_add_pd(
            xi_209, _mm256_set_pd(
                        ((points_down && ctr_1 <= 0)
                             ? (1.0)
                             : ((points_up && ctr_1 >= 63) ? (-1.0) : (0.0))),
                        ((points_down && ctr_1 <= 0)
                             ? (1.0)
                             : ((points_up && ctr_1 >= 63) ? (-1.0) : (0.0))),
                        ((points_down && ctr_1 <= 0)
                             ? (1.0)
                             : ((points_up && ctr_1 >= 63) ? (-1.0) : (0.0))),
                        ((points_down && ctr_1 <= 0)
                             ? (1.0)
                             : ((points_up && ctr_1 >= 63) ? (-1.0) : (0.0)))));
        const __m256d xi_0 = _mm256_mul_pd(u_0, xi_210);
        const __m256d xi_1 = _mm256_mul_pd(
            xi_0, _mm256_set_pd(0.333333333333333, 0.333333333333333,
                                0.333333333333333, 0.333333333333333));
        const __m256d xi_7 =
            _mm256_mul_pd(xi_1, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_65 = _mm256_mul_pd(rho, (_mm256_mul_pd(u_0, u_0)));
        const __m256d xi_128 = _mm256_mul_pd(rho, u_0);
        const __m256d xi_131 = _mm256_add_pd(
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(_mm256_add_pd(xi_128, xi_130), xi_211),
                    xi_217),
                xi_222),
            xi_93);
        const __m256d xi_132 = _mm256_mul_pd(
            xi_131, _mm256_set_pd(xi_116, xi_116, xi_116, xi_116));
        const __m256d xi_159 = _mm256_mul_pd(
            xi_131, _mm256_set_pd(xi_156, xi_156, xi_156, xi_156));
        const __m256d u_1 = xi_224;
        const __m256d xi_2 = _mm256_mul_pd(u_1, xi_215);
        const __m256d xi_3 = _mm256_mul_pd(
            xi_2, _mm256_set_pd(0.333333333333333, 0.333333333333333,
                                0.333333333333333, 0.333333333333333));
        const __m256d xi_8 =
            _mm256_mul_pd(xi_3, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_29 =
            _mm256_mul_pd(u_1, _mm256_set_pd(0.5, 0.5, 0.5, 0.5));
        const __m256d xi_32 =
            _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(u_0, xi_30),
                                        _mm256_mul_pd(xi_210, xi_29)),
                          _mm256_set_pd(xi_31, xi_31, xi_31, xi_31));
        const __m256d xi_33 =
            _mm256_mul_pd(xi_32, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_68 = _mm256_mul_pd(rho, (_mm256_mul_pd(u_1, u_1)));
        const __m256d xi_71 = _mm256_add_pd(_mm256_add_pd(xi_67, xi_68), xi_70);
        const __m256d xi_113 = _mm256_mul_pd(rho, u_1);
        const __m256d xi_115 = _mm256_add_pd(
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(_mm256_add_pd(xi_101, xi_113), xi_114),
                    xi_218),
                xi_222),
            xi_96);
        const __m256d xi_117 = _mm256_mul_pd(
            xi_115, _mm256_set_pd(xi_116, xi_116, xi_116, xi_116));
        const __m256d xi_152 = _mm256_mul_pd(
            _mm256_add_pd(
                _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(u_0, xi_113), xi_114),
                              xi_211),
                xi_67),
            _mm256_set_pd(xi_151, xi_151, xi_151, xi_151));
        const __m256d xi_157 = _mm256_mul_pd(
            xi_115, _mm256_set_pd(xi_156, xi_156, xi_156, xi_156));
        const __m256d xi_158 = _mm256_add_pd(xi_111, xi_157);
        const __m256d xi_167 =
            _mm256_mul_pd(xi_157, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_168 = _mm256_add_pd(xi_112, xi_167);
        const __m256d xi_183 = _mm256_add_pd(
            _mm256_add_pd(
                _mm256_mul_pd(xi_181, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                xi_157),
            xi_182);
        const __m256d xi_184 = _mm256_add_pd(
            _mm256_add_pd(
                _mm256_mul_pd(xi_182, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                xi_167),
            xi_181);
        const __m256d u_2 = xi_206;
        const __m256d xi_4 = _mm256_mul_pd(u_2, xi_219);
        const __m256d xi_5 = _mm256_mul_pd(
            xi_4, _mm256_set_pd(0.333333333333333, 0.333333333333333,
                                0.333333333333333, 0.333333333333333));
        const __m256d xi_6 =
            _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(xi_1, xi_3), xi_5),
                          _mm256_set_pd(-omega_bulk + 2.0, -omega_bulk + 2.0,
                                        -omega_bulk + 2.0, -omega_bulk + 2.0));
        const __m256d xi_9 = _mm256_add_pd(
            _mm256_add_pd(_mm256_mul_pd(xi_4, _mm256_set_pd(0.666666666666667,
                                                            0.666666666666667,
                                                            0.666666666666667,
                                                            0.666666666666667)),
                          xi_7),
            xi_8);
        const __m256d xi_12 =
            _mm256_mul_pd(xi_5, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_13 = _mm256_add_pd(
            _mm256_add_pd(_mm256_mul_pd(xi_2, _mm256_set_pd(0.666666666666667,
                                                            0.666666666666667,
                                                            0.666666666666667,
                                                            0.666666666666667)),
                          xi_12),
            xi_7);
        const __m256d xi_14 = _mm256_add_pd(
            _mm256_add_pd(_mm256_mul_pd(xi_0, _mm256_set_pd(0.666666666666667,
                                                            0.666666666666667,
                                                            0.666666666666667,
                                                            0.666666666666667)),
                          xi_12),
            xi_8);
        const __m256d xi_17 =
            _mm256_mul_pd(xi_9, _mm256_set_pd(xi_16, xi_16, xi_16, xi_16));
        const __m256d xi_18 =
            _mm256_mul_pd(xi_17, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_19 =
            _mm256_mul_pd(xi_14, _mm256_set_pd(xi_16, xi_16, xi_16, xi_16));
        const __m256d xi_20 =
            _mm256_mul_pd(xi_19, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_22 = _mm256_add_pd(
            _mm256_add_pd(
                _mm256_mul_pd(xi_13, _mm256_set_pd(xi_21, xi_21, xi_21, xi_21)),
                xi_18),
            xi_20);
        const __m256d xi_24 =
            _mm256_mul_pd(xi_13, _mm256_set_pd(xi_16, xi_16, xi_16, xi_16));
        const __m256d xi_25 =
            _mm256_mul_pd(xi_24, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_26 = _mm256_add_pd(
            _mm256_add_pd(
                _mm256_mul_pd(xi_14, _mm256_set_pd(xi_21, xi_21, xi_21, xi_21)),
                xi_18),
            xi_25);
        const __m256d xi_28 = _mm256_add_pd(
            _mm256_add_pd(
                _mm256_mul_pd(xi_9, _mm256_set_pd(xi_21, xi_21, xi_21, xi_21)),
                xi_20),
            xi_25);
        const __m256d xi_35 = _mm256_add_pd(
            _mm256_mul_pd(xi_34, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_19);
        const __m256d xi_37 = _mm256_mul_pd(
            _mm256_mul_pd(xi_9, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
            _mm256_set_pd(xi_36, xi_36, xi_36, xi_36));
        const __m256d xi_39 =
            _mm256_mul_pd(xi_6, _mm256_set_pd(0.125, 0.125, 0.125, 0.125));
        const __m256d xi_40 = _mm256_add_pd(xi_24, xi_39);
        const __m256d xi_41 = _mm256_add_pd(xi_38, xi_40);
        const __m256d xi_42 = _mm256_add_pd(xi_37, xi_41);
        const __m256d xi_43 = _mm256_add_pd(xi_19, xi_34);
        const __m256d xi_44 = _mm256_add_pd(
            _mm256_mul_pd(xi_38, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_40);
        const __m256d xi_45 = _mm256_add_pd(xi_37, xi_44);
        const __m256d xi_46 =
            _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(u_2, xi_30),
                                        _mm256_mul_pd(xi_219, xi_29)),
                          _mm256_set_pd(xi_31, xi_31, xi_31, xi_31));
        const __m256d xi_47 = _mm256_mul_pd(
            _mm256_mul_pd(xi_14, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
            _mm256_set_pd(xi_36, xi_36, xi_36, xi_36));
        const __m256d xi_49 = _mm256_add_pd(xi_17, xi_48);
        const __m256d xi_50 = _mm256_add_pd(xi_47, xi_49);
        const __m256d xi_51 =
            _mm256_mul_pd(xi_46, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_52 = _mm256_mul_pd(
            _mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(u_0, xi_219),
                                        _mm256_set_pd(0.5, 0.5, 0.5, 0.5)),
                          _mm256_mul_pd(_mm256_mul_pd(u_2, xi_210),
                                        _mm256_set_pd(0.5, 0.5, 0.5, 0.5))),
            _mm256_set_pd(xi_31, xi_31, xi_31, xi_31));
        const __m256d xi_53 =
            _mm256_mul_pd(xi_52, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_54 = _mm256_mul_pd(
            _mm256_mul_pd(xi_13, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
            _mm256_set_pd(xi_36, xi_36, xi_36, xi_36));
        const __m256d xi_55 = _mm256_add_pd(_mm256_add_pd(xi_39, xi_49), xi_54);
        const __m256d xi_56 = _mm256_add_pd(
            _mm256_mul_pd(xi_48, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_17);
        const __m256d xi_57 = _mm256_add_pd(xi_47, xi_56);
        const __m256d xi_58 = _mm256_add_pd(_mm256_add_pd(xi_39, xi_54), xi_56);
        const __m256d xi_66 = _mm256_mul_pd(rho, (_mm256_mul_pd(u_2, u_2)));
        const __m256d xi_86 = _mm256_mul_pd(
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_add_pd(_mm256_add_pd(xi_204, xi_65), xi_66),
                        xi_71),
                    xi_80),
                xi_85),
            _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk, omega_bulk));
        const __m256d xi_119 = _mm256_add_pd(
            _mm256_mul_pd(xi_66, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
            xi_202);
        const __m256d xi_120 = _mm256_mul_pd(
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_add_pd(_mm256_add_pd(xi_118, xi_119), xi_200),
                        xi_223),
                    xi_71),
                xi_97),
            _mm256_set_pd(omega_shear, omega_shear, omega_shear, omega_shear));
        const __m256d xi_121 =
            _mm256_mul_pd(xi_120, _mm256_set_pd(0.125, 0.125, 0.125, 0.125));
        const __m256d xi_124 = _mm256_mul_pd(
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
                                                            _mm256_add_pd(
                                                                _mm256_mul_pd(
                                                                    xi_217,
                                                                    _mm256_set_pd(
                                                                        -2.0,
                                                                        -2.0,
                                                                        -2.0,
                                                                        -2.0)),
                                                                _mm256_mul_pd(
                                                                    xi_221,
                                                                    _mm256_set_pd(
                                                                        -2.0,
                                                                        -2.0,
                                                                        -2.0,
                                                                        -2.0))),
                                                            _mm256_mul_pd(
                                                                xi_65,
                                                                _mm256_set_pd(
                                                                    2.0, 2.0,
                                                                    2.0, 2.0))),
                                                        _mm256_mul_pd(
                                                            xi_68,
                                                            _mm256_set_pd(
                                                                -1.0, -1.0,
                                                                -1.0, -1.0))),
                                                    xi_119),
                                                xi_201),
                                            xi_218),
                                        xi_220),
                                    xi_62),
                                xi_67),
                            xi_70),
                        xi_72),
                    xi_73),
                xi_85),
            _mm256_set_pd(omega_shear, omega_shear, omega_shear, omega_shear));
        const __m256d xi_125 = _mm256_mul_pd(
            xi_124, _mm256_set_pd(-0.0416666666666667, -0.0416666666666667,
                                  -0.0416666666666667, -0.0416666666666667));
        const __m256d xi_126 = _mm256_add_pd(
            _mm256_mul_pd(xi_61, _mm256_set_pd(-0.05, -0.05, -0.05, -0.05)),
            xi_125);
        const __m256d xi_127 = _mm256_add_pd(
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_mul_pd(xi_64, _mm256_set_pd(0.0142857142857143,
                                                       0.0142857142857143,
                                                       0.0142857142857143,
                                                       0.0142857142857143)),
                    xi_121),
                xi_123),
            xi_126);
        const __m256d xi_139 = _mm256_add_pd(
            _mm256_add_pd(
                _mm256_mul_pd(xi_124, _mm256_set_pd(0.0833333333333333,
                                                    0.0833333333333333,
                                                    0.0833333333333333,
                                                    0.0833333333333333)),
                _mm256_mul_pd(xi_64, _mm256_set_pd(-0.0357142857142857,
                                                   -0.0357142857142857,
                                                   -0.0357142857142857,
                                                   -0.0357142857142857))),
            xi_123);
        const __m256d xi_148 = _mm256_add_pd(
            _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(rho, u_2), xi_118),
                          xi_142),
            xi_94);
        const __m256d xi_149 = _mm256_mul_pd(
            xi_148, _mm256_set_pd(xi_116, xi_116, xi_116, xi_116));
        const __m256d xi_150 = _mm256_add_pd(
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_mul_pd(xi_121,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_pd(xi_64,
                                      _mm256_set_pd(-0.0214285714285714,
                                                    -0.0214285714285714,
                                                    -0.0214285714285714,
                                                    -0.0214285714285714))),
                    _mm256_mul_pd(xi_90, _mm256_set_pd(0.0158730158730159,
                                                       0.0158730158730159,
                                                       0.0158730158730159,
                                                       0.0158730158730159))),
                xi_122),
            xi_126);
        const __m256d xi_153 = _mm256_mul_pd(
            xi_120, _mm256_set_pd(0.0625, 0.0625, 0.0625, 0.0625));
        const __m256d xi_155 = _mm256_add_pd(
            _mm256_add_pd(
                _mm256_mul_pd(xi_152, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                xi_153),
            xi_154);
        const __m256d xi_160 = _mm256_add_pd(
            _mm256_mul_pd(
                xi_86, _mm256_set_pd(0.0416666666666667, 0.0416666666666667,
                                     0.0416666666666667, 0.0416666666666667)),
            _mm256_mul_pd(
                xi_98, _mm256_set_pd(0.0277777777777778, 0.0277777777777778,
                                     0.0277777777777778, 0.0277777777777778)));
        const __m256d xi_161 = _mm256_add_pd(
            _mm256_mul_pd(
                xi_124, _mm256_set_pd(0.0208333333333333, 0.0208333333333333,
                                      0.0208333333333333, 0.0208333333333333)),
            xi_160);
        const __m256d xi_162 = _mm256_add_pd(
            _mm256_mul_pd(xi_159, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
            xi_161);
        const __m256d xi_163 = _mm256_add_pd(xi_140, xi_162);
        const __m256d xi_164 =
            _mm256_add_pd(_mm256_add_pd(xi_152, xi_153), xi_154);
        const __m256d xi_165 = _mm256_add_pd(xi_159, xi_161);
        const __m256d xi_166 = _mm256_add_pd(xi_138, xi_165);
        const __m256d xi_170 = _mm256_mul_pd(
            _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(u_2, xi_113), xi_100),
                          xi_107),
            _mm256_set_pd(xi_151, xi_151, xi_151, xi_151));
        const __m256d xi_172 = _mm256_add_pd(
            _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_125, xi_160), xi_169),
                          xi_170),
            xi_171);
        const __m256d xi_178 = _mm256_mul_pd(
            xi_148, _mm256_set_pd(xi_156, xi_156, xi_156, xi_156));
        const __m256d xi_179 = _mm256_add_pd(xi_177, xi_178);
        const __m256d xi_180 = _mm256_add_pd(
            _mm256_add_pd(
                _mm256_mul_pd(xi_174, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                xi_176),
            xi_179);
        const __m256d xi_185 = _mm256_add_pd(
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_mul_pd(xi_170,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                        xi_125),
                    xi_160),
                xi_169),
            xi_171);
        const __m256d xi_186 = _mm256_mul_pd(
            _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(u_2, xi_128), xi_129),
                          xi_136),
            _mm256_set_pd(xi_151, xi_151, xi_151, xi_151));
        const __m256d xi_187 =
            _mm256_mul_pd(xi_153, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
        const __m256d xi_189 = _mm256_add_pd(
            _mm256_add_pd(
                _mm256_mul_pd(xi_186, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                xi_187),
            xi_188);
        const __m256d xi_190 = _mm256_add_pd(xi_146, xi_179);
        const __m256d xi_193 = _mm256_add_pd(
            _mm256_add_pd(
                _mm256_mul_pd(xi_191, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                xi_162),
            xi_192);
        const __m256d xi_194 =
            _mm256_add_pd(_mm256_add_pd(xi_186, xi_187), xi_188);
        const __m256d xi_195 = _mm256_add_pd(
            _mm256_add_pd(
                _mm256_mul_pd(xi_192, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                xi_165),
            xi_191);
        const __m256d xi_196 = _mm256_add_pd(
            _mm256_mul_pd(xi_178, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
            xi_177);
        const __m256d xi_197 = _mm256_add_pd(
            _mm256_add_pd(
                _mm256_mul_pd(xi_176, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                xi_174),
            xi_196);
        const __m256d xi_198 = _mm256_add_pd(xi_147, xi_196);
        const __m256d forceTerm_0 = _mm256_add_pd(
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_mul_pd(xi_6, _mm256_set_pd(-1.5, -1.5, -1.5, -1.5)),
                    _mm256_mul_pd(
                        _mm256_mul_pd(xi_13,
                                      _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                        _mm256_set_pd(xi_11, xi_11, xi_11, xi_11))),
                _mm256_mul_pd(
                    _mm256_mul_pd(xi_14, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                    _mm256_set_pd(xi_11, xi_11, xi_11, xi_11))),
            _mm256_mul_pd(
                _mm256_mul_pd(xi_9, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                _mm256_set_pd(xi_11, xi_11, xi_11, xi_11)));
        const __m256d forceTerm_1 = _mm256_add_pd(xi_15, xi_22);
        const __m256d forceTerm_2 = _mm256_add_pd(
            _mm256_mul_pd(xi_15, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_22);
        const __m256d forceTerm_3 = _mm256_add_pd(
            _mm256_mul_pd(xi_23, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_26);
        const __m256d forceTerm_4 = _mm256_add_pd(xi_23, xi_26);
        const __m256d forceTerm_5 = _mm256_add_pd(xi_27, xi_28);
        const __m256d forceTerm_6 = _mm256_add_pd(
            _mm256_mul_pd(xi_27, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_28);
        const __m256d forceTerm_7 =
            _mm256_add_pd(_mm256_add_pd(xi_33, xi_35), xi_42);
        const __m256d forceTerm_8 =
            _mm256_add_pd(_mm256_add_pd(xi_32, xi_42), xi_43);
        const __m256d forceTerm_9 =
            _mm256_add_pd(_mm256_add_pd(xi_32, xi_35), xi_45);
        const __m256d forceTerm_10 =
            _mm256_add_pd(_mm256_add_pd(xi_33, xi_43), xi_45);
        const __m256d forceTerm_11 =
            _mm256_add_pd(_mm256_add_pd(xi_41, xi_46), xi_50);
        const __m256d forceTerm_12 =
            _mm256_add_pd(_mm256_add_pd(xi_44, xi_50), xi_51);
        const __m256d forceTerm_13 =
            _mm256_add_pd(_mm256_add_pd(xi_35, xi_53), xi_55);
        const __m256d forceTerm_14 =
            _mm256_add_pd(_mm256_add_pd(xi_43, xi_52), xi_55);
        const __m256d forceTerm_15 =
            _mm256_add_pd(_mm256_add_pd(xi_41, xi_51), xi_57);
        const __m256d forceTerm_16 =
            _mm256_add_pd(_mm256_add_pd(xi_44, xi_46), xi_57);
        const __m256d forceTerm_17 =
            _mm256_add_pd(_mm256_add_pd(xi_35, xi_52), xi_58);
        const __m256d forceTerm_18 =
            _mm256_add_pd(_mm256_add_pd(xi_43, xi_53), xi_58);
        _mm256_store_pd(
            &_data_pdfs_20_30_10[ctr_0],
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_add_pd(
                            _mm256_add_pd(
                                _mm256_add_pd(
                                    _mm256_mul_pd(
                                        xi_61,
                                        _mm256_set_pd(0.1, 0.1, 0.1, 0.1)),
                                    _mm256_mul_pd(
                                        xi_64,
                                        _mm256_set_pd(0.0428571428571429,
                                                      0.0428571428571429,
                                                      0.0428571428571429,
                                                      0.0428571428571429))),
                                _mm256_mul_pd(
                                    xi_86,
                                    _mm256_set_pd(-0.5, -0.5, -0.5, -0.5))),
                            _mm256_mul_pd(xi_90,
                                          _mm256_set_pd(0.0238095238095238,
                                                        0.0238095238095238,
                                                        0.0238095238095238,
                                                        0.0238095238095238))),
                        _mm256_mul_pd(xi_98, _mm256_set_pd(0.333333333333333,
                                                           0.333333333333333,
                                                           0.333333333333333,
                                                           0.333333333333333))),
                    forceTerm_0),
                xi_204));
        _mm256_store_pd(
            &_data_pdfs_20_31_10[ctr_0],
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_add_pd(
                            _mm256_add_pd(
                                _mm256_mul_pd(
                                    xi_104,
                                    _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                                forceTerm_1),
                            xi_112),
                        xi_117),
                    xi_127),
                xi_220));
        _mm256_store_pd(
            &_data_pdfs_20_32_10[ctr_0],
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_add_pd(
                            _mm256_add_pd(
                                _mm256_mul_pd(
                                    xi_117,
                                    _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                                forceTerm_2),
                            xi_104),
                        xi_111),
                    xi_127),
                xi_218));
        _mm256_store_pd(
            &_data_pdfs_20_33_10[ctr_0],
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_add_pd(
                            _mm256_add_pd(
                                _mm256_mul_pd(
                                    xi_132,
                                    _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                                forceTerm_3),
                            xi_135),
                        xi_138),
                    xi_139),
                xi_217));
        _mm256_store_pd(
            &_data_pdfs_20_34_10[ctr_0],
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_add_pd(
                            _mm256_add_pd(
                                _mm256_mul_pd(
                                    xi_135,
                                    _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                                forceTerm_4),
                            xi_132),
                        xi_139),
                    xi_140),
                xi_221));
        _mm256_store_pd(
            &_data_pdfs_20_35_10[ctr_0],
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_add_pd(
                            _mm256_add_pd(
                                _mm256_mul_pd(
                                    xi_144,
                                    _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                                forceTerm_5),
                            xi_147),
                        xi_149),
                    xi_150),
                xi_202));
        _mm256_store_pd(
            &_data_pdfs_20_36_10[ctr_0],
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(
                        _mm256_add_pd(
                            _mm256_add_pd(
                                _mm256_mul_pd(
                                    xi_149,
                                    _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)),
                                forceTerm_6),
                            xi_144),
                        xi_146),
                    xi_150),
                xi_201));
        _mm256_store_pd(
            &_data_pdfs_20_37_10[ctr_0],
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(_mm256_add_pd(forceTerm_7, xi_155), xi_158),
                    xi_163),
                xi_211));
        _mm256_store_pd(
            &_data_pdfs_20_38_10[ctr_0],
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(_mm256_add_pd(forceTerm_8, xi_158), xi_164),
                    xi_166),
                xi_208));
        _mm256_store_pd(
            &_data_pdfs_20_39_10[ctr_0],
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(_mm256_add_pd(forceTerm_9, xi_163), xi_164),
                    xi_168),
                xi_222));
        _mm256_store_pd(
            &_data_pdfs_20_310_10[ctr_0],
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(_mm256_add_pd(forceTerm_10, xi_155), xi_166),
                    xi_168),
                xi_203));
        _mm256_store_pd(
            &_data_pdfs_20_311_10[ctr_0],
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(_mm256_add_pd(forceTerm_11, xi_172), xi_180),
                    xi_183),
                xi_214));
        _mm256_store_pd(
            &_data_pdfs_20_312_10[ctr_0],
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(_mm256_add_pd(forceTerm_12, xi_180), xi_184),
                    xi_185),
                xi_199));
        _mm256_store_pd(
            &_data_pdfs_20_313_10[ctr_0],
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(_mm256_add_pd(forceTerm_13, xi_189), xi_190),
                    xi_193),
                xi_200));
        _mm256_store_pd(
            &_data_pdfs_20_314_10[ctr_0],
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(_mm256_add_pd(forceTerm_14, xi_190), xi_194),
                    xi_195),
                xi_223));
        _mm256_store_pd(
            &_data_pdfs_20_315_10[ctr_0],
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(_mm256_add_pd(forceTerm_15, xi_183), xi_185),
                    xi_197),
                xi_213));
        _mm256_store_pd(
            &_data_pdfs_20_316_10[ctr_0],
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(_mm256_add_pd(forceTerm_16, xi_172), xi_184),
                    xi_197),
                xi_207));
        _mm256_store_pd(
            &_data_pdfs_20_317_10[ctr_0],
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(_mm256_add_pd(forceTerm_17, xi_193), xi_194),
                    xi_198),
                xi_216));
        _mm256_store_pd(
            &_data_pdfs_20_318_10[ctr_0],
            _mm256_add_pd(
                _mm256_add_pd(
                    _mm256_add_pd(_mm256_add_pd(forceTerm_18, xi_189), xi_195),
                    xi_198),
                xi_212));
      }
    }
  }
}
} // namespace internal_collidesweepleesedwardsavx

void CollideSweepLeesEdwardsAVX::operator()(IBlock *block) {
  auto force = block->getData<field::GhostLayerField<double, 3>>(forceID);
  auto density = block->getData<field::GhostLayerField<double, 1>>(densityID);
  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);
  auto velocity = block->getData<field::GhostLayerField<double, 3>>(velocityID);

  auto &points_up = this->points_up_;
  auto &omega_shear = this->omega_shear_;
  auto &omega_bulk = this->omega_bulk_;
  auto &omega_even = this->omega_even_;
  auto &omega_odd = this->omega_odd_;
  auto &points_down = this->points_down_;
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(density->nrOfGhostLayers()));
  double *RESTRICT const _data_density = density->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(density->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force->nrOfGhostLayers()));
  double *RESTRICT const _data_force = force->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()));
  double *RESTRICT _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(velocity->nrOfGhostLayers()));
  double *RESTRICT const _data_velocity = velocity->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(density->xSizeWithGhostLayer(),
                                int64_t(cell_idx_c(density->xSize()) + 0));
  const int64_t _size_density_0 = int64_t(cell_idx_c(density->xSize()) + 0);
  WALBERLA_ASSERT_EQUAL(density->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(density->ySizeWithGhostLayer(),
                                int64_t(cell_idx_c(density->ySize()) + 0));
  const int64_t _size_density_1 = int64_t(cell_idx_c(density->ySize()) + 0);
  WALBERLA_ASSERT_EQUAL(density->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(density->zSizeWithGhostLayer(),
                                int64_t(cell_idx_c(density->zSize()) + 0));
  const int64_t _size_density_2 = int64_t(cell_idx_c(density->zSize()) + 0);
  WALBERLA_ASSERT_EQUAL(density->layout(), field::fzyx);
  const int64_t _stride_density_1 = int64_t(density->yStride());
  const int64_t _stride_density_2 = int64_t(density->zStride());
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
  const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
  const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
  internal_collidesweepleesedwardsavx::collidesweepleesedwardsavx(
      _data_density, _data_force, _data_pdfs, _data_velocity, _size_density_0,
      _size_density_1, _size_density_2, _stride_density_1, _stride_density_2,
      _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_1,
      _stride_pdfs_2, _stride_pdfs_3, _stride_velocity_1, _stride_velocity_2,
      _stride_velocity_3, omega_bulk, omega_even, omega_odd, omega_shear,
      points_down, points_up);
}

void CollideSweepLeesEdwardsAVX::runOnCellInterval(
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
  auto density = block->getData<field::GhostLayerField<double, 1>>(densityID);
  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);
  auto velocity = block->getData<field::GhostLayerField<double, 3>>(velocityID);

  auto &points_up = this->points_up_;
  auto &omega_shear = this->omega_shear_;
  auto &omega_bulk = this->omega_bulk_;
  auto &omega_even = this->omega_even_;
  auto &omega_odd = this->omega_odd_;
  auto &points_down = this->points_down_;
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(density->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(density->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(density->nrOfGhostLayers()));
  double *RESTRICT const _data_density =
      density->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(density->layout(), field::fzyx);
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
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(velocity->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(velocity->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(velocity->nrOfGhostLayers()));
  double *RESTRICT const _data_velocity =
      velocity->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(density->xSizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.xSize()) + 0));
  const int64_t _size_density_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
  WALBERLA_ASSERT_EQUAL(density->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(density->ySizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.ySize()) + 0));
  const int64_t _size_density_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
  WALBERLA_ASSERT_EQUAL(density->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(density->zSizeWithGhostLayer(),
                                int64_t(cell_idx_c(ci.zSize()) + 0));
  const int64_t _size_density_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
  WALBERLA_ASSERT_EQUAL(density->layout(), field::fzyx);
  const int64_t _stride_density_1 = int64_t(density->yStride());
  const int64_t _stride_density_2 = int64_t(density->zStride());
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
  const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
  const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
  internal_collidesweepleesedwardsavx::collidesweepleesedwardsavx(
      _data_density, _data_force, _data_pdfs, _data_velocity, _size_density_0,
      _size_density_1, _size_density_2, _stride_density_1, _stride_density_2,
      _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_1,
      _stride_pdfs_2, _stride_pdfs_3, _stride_velocity_1, _stride_velocity_2,
      _stride_velocity_3, omega_bulk, omega_even, omega_odd, omega_shear,
      points_down, points_up);
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