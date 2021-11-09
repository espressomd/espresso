// kernel generated with pystencils v0.3.4+4.g4fecf0c, lbmpy v0.3.4+6.g2faceda,
// lbmpy_walberla/pystencils_walberla from commit
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
//! \\file CollideSweepSinglePrecisionThermalizedAVX.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

#include <cmath>

#include "CollideSweepSinglePrecisionThermalizedAVX.h"
#include "core/DataTypes.h"
#include "core/Macros.h"

#include "philox_rand.h"

#include <immintrin.h>

#define FUNC_PREFIX

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) ||                                  \
    (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning push
#pragma warning(disable : 1599)
#endif

using namespace std;

namespace walberla {
namespace pystencils {

namespace internal_collidesweepsingleprecisionthermalizedavx {
static FUNC_PREFIX void collidesweepsingleprecisionthermalizedavx(
    float *RESTRICT const _data_force, float *RESTRICT _data_pdfs,
    int64_t const _size_force_0, int64_t const _size_force_1,
    int64_t const _size_force_2, int64_t const _stride_force_1,
    int64_t const _stride_force_2, int64_t const _stride_force_3,
    int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
    int64_t const _stride_pdfs_3, uint32_t block_offset_0,
    uint32_t block_offset_1, uint32_t block_offset_2, double kT,
    double omega_bulk, double omega_even, double omega_odd, double omega_shear,
    uint32_t seed, uint32_t time_step) {
  const double xi_25 = -omega_bulk;
  const double xi_36 = -omega_shear;
  const double xi_37 = xi_36 + 2.0;
  const double xi_38 = xi_37 * 0.5;
  const double xi_43 = xi_37 * 0.0833333333333333;
  const double xi_48 = xi_37 * 0.166666666666667;
  const double xi_58 = xi_37 * 0.25;
  const double xi_63 = xi_37 * 0.0416666666666667;
  const float xi_90 = 2.4494897427831779;
  const double xi_115 = omega_odd * 0.25;
  const double xi_131 = omega_odd * 0.0833333333333333;
  const double xi_196 = omega_shear * 0.25;
  const double xi_211 = omega_odd * 0.0416666666666667;
  const double xi_213 = omega_odd * 0.125;
  const int64_t rr_0 = 0.0;
  const float xi_120 = rr_0 * 0.166666666666667;
  const float xi_186 = rr_0 * 0.0833333333333333;
  for (int64_t ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1) {
    float *RESTRICT _data_force_20_31 =
        _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    float *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_312 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_315 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_318 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_35 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_32 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_313 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_34 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    float *RESTRICT _data_pdfs_20_317 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_32 =
        _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    float *RESTRICT _data_pdfs_20_31 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_33 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_36 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_311 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_314 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_316 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    for (int64_t ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1) {
      float *RESTRICT _data_force_20_31_10 =
          _stride_force_1 * ctr_1 + _data_force_20_31;
      float *RESTRICT _data_pdfs_20_37_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_37;
      float *RESTRICT _data_pdfs_20_312_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_312;
      float *RESTRICT _data_pdfs_20_315_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_315;
      float *RESTRICT _data_pdfs_20_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_318;
      float *RESTRICT _data_pdfs_20_35_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_35;
      float *RESTRICT _data_pdfs_20_32_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_32;
      float *RESTRICT _data_pdfs_20_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_313;
      float *RESTRICT _data_pdfs_20_310_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_310;
      float *RESTRICT _data_pdfs_20_34_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      float *RESTRICT _data_pdfs_20_38_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_38;
      float *RESTRICT _data_pdfs_20_30_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      float *RESTRICT _data_pdfs_20_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_317;
      float *RESTRICT _data_force_20_32_10 =
          _stride_force_1 * ctr_1 + _data_force_20_32;
      float *RESTRICT _data_pdfs_20_31_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_31;
      float *RESTRICT _data_pdfs_20_33_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      float *RESTRICT _data_pdfs_20_36_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_36;
      float *RESTRICT _data_pdfs_20_311_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_311;
      float *RESTRICT _data_pdfs_20_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_314;
      float *RESTRICT _data_pdfs_20_316_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_316;
      float *RESTRICT _data_pdfs_20_39_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_39;
      float *RESTRICT _data_force_20_30_10 =
          _stride_force_1 * ctr_1 + _data_force_20_30;
      for (int64_t ctr_0 = 0;
           ctr_0 < ((_size_force_0) % (8) == 0
                        ? _size_force_0
                        : ((int64_t)((_size_force_0) / (8)) + 1) * (8));
           ctr_0 += 8) {
        const __m256 xi_248 = _mm256_load_ps(&_data_force_20_31_10[ctr_0]);
        const __m256 xi_249 = _mm256_load_ps(&_data_pdfs_20_37_10[ctr_0]);
        const __m256 xi_250 = _mm256_load_ps(&_data_pdfs_20_312_10[ctr_0]);
        const __m256 xi_251 = _mm256_load_ps(&_data_pdfs_20_315_10[ctr_0]);
        const __m256 xi_252 = _mm256_load_ps(&_data_pdfs_20_318_10[ctr_0]);
        const __m256 xi_253 = _mm256_load_ps(&_data_pdfs_20_35_10[ctr_0]);
        const __m256 xi_254 = _mm256_load_ps(&_data_pdfs_20_32_10[ctr_0]);
        const __m256 xi_255 = _mm256_load_ps(&_data_pdfs_20_313_10[ctr_0]);
        const __m256 xi_256 = _mm256_load_ps(&_data_pdfs_20_310_10[ctr_0]);
        const __m256 xi_257 = _mm256_load_ps(&_data_pdfs_20_34_10[ctr_0]);
        const __m256 xi_258 = _mm256_load_ps(&_data_pdfs_20_38_10[ctr_0]);
        const __m256 xi_259 = _mm256_load_ps(&_data_pdfs_20_30_10[ctr_0]);
        const __m256 xi_260 = _mm256_load_ps(&_data_pdfs_20_317_10[ctr_0]);
        const __m256 xi_261 = _mm256_load_ps(&_data_force_20_32_10[ctr_0]);
        const __m256 xi_262 = _mm256_load_ps(&_data_pdfs_20_31_10[ctr_0]);
        const __m256 xi_263 = _mm256_load_ps(&_data_pdfs_20_33_10[ctr_0]);
        const __m256 xi_264 = _mm256_load_ps(&_data_pdfs_20_36_10[ctr_0]);
        const __m256 xi_265 = _mm256_load_ps(&_data_pdfs_20_311_10[ctr_0]);
        const __m256 xi_266 = _mm256_load_ps(&_data_pdfs_20_314_10[ctr_0]);
        const __m256 xi_267 = _mm256_load_ps(&_data_pdfs_20_316_10[ctr_0]);
        const __m256 xi_268 = _mm256_load_ps(&_data_pdfs_20_39_10[ctr_0]);
        const __m256 xi_269 = _mm256_load_ps(&_data_force_20_30_10[ctr_0]);

        __m256d random_7_0;
        __m256d random_7_1;
        philox_double2(time_step,
                       _mm256_add_epi32(
                           _mm256_add_epi32(
                               _mm256_set_epi32(block_offset_0, block_offset_0,
                                                block_offset_0, block_offset_0,
                                                block_offset_0, block_offset_0,
                                                block_offset_0, block_offset_0),
                               _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0)),
                           _mm256_set_epi32(ctr_0, ctr_0, ctr_0, ctr_0, ctr_0,
                                            ctr_0, ctr_0, ctr_0)),
                       block_offset_1 + ctr_1, block_offset_2 + ctr_2, 7, seed,
                       random_7_0, random_7_1);

        __m256d random_6_0;
        __m256d random_6_1;
        philox_double2(time_step,
                       _mm256_add_epi32(
                           _mm256_add_epi32(
                               _mm256_set_epi32(block_offset_0, block_offset_0,
                                                block_offset_0, block_offset_0,
                                                block_offset_0, block_offset_0,
                                                block_offset_0, block_offset_0),
                               _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0)),
                           _mm256_set_epi32(ctr_0, ctr_0, ctr_0, ctr_0, ctr_0,
                                            ctr_0, ctr_0, ctr_0)),
                       block_offset_1 + ctr_1, block_offset_2 + ctr_2, 6, seed,
                       random_6_0, random_6_1);

        __m256d random_5_0;
        __m256d random_5_1;
        philox_double2(time_step,
                       _mm256_add_epi32(
                           _mm256_add_epi32(
                               _mm256_set_epi32(block_offset_0, block_offset_0,
                                                block_offset_0, block_offset_0,
                                                block_offset_0, block_offset_0,
                                                block_offset_0, block_offset_0),
                               _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0)),
                           _mm256_set_epi32(ctr_0, ctr_0, ctr_0, ctr_0, ctr_0,
                                            ctr_0, ctr_0, ctr_0)),
                       block_offset_1 + ctr_1, block_offset_2 + ctr_2, 5, seed,
                       random_5_0, random_5_1);

        __m256d random_4_0;
        __m256d random_4_1;
        philox_double2(time_step,
                       _mm256_add_epi32(
                           _mm256_add_epi32(
                               _mm256_set_epi32(block_offset_0, block_offset_0,
                                                block_offset_0, block_offset_0,
                                                block_offset_0, block_offset_0,
                                                block_offset_0, block_offset_0),
                               _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0)),
                           _mm256_set_epi32(ctr_0, ctr_0, ctr_0, ctr_0, ctr_0,
                                            ctr_0, ctr_0, ctr_0)),
                       block_offset_1 + ctr_1, block_offset_2 + ctr_2, 4, seed,
                       random_4_0, random_4_1);

        __m256d random_3_0;
        __m256d random_3_1;
        philox_double2(time_step,
                       _mm256_add_epi32(
                           _mm256_add_epi32(
                               _mm256_set_epi32(block_offset_0, block_offset_0,
                                                block_offset_0, block_offset_0,
                                                block_offset_0, block_offset_0,
                                                block_offset_0, block_offset_0),
                               _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0)),
                           _mm256_set_epi32(ctr_0, ctr_0, ctr_0, ctr_0, ctr_0,
                                            ctr_0, ctr_0, ctr_0)),
                       block_offset_1 + ctr_1, block_offset_2 + ctr_2, 3, seed,
                       random_3_0, random_3_1);

        __m256d random_2_0;
        __m256d random_2_1;
        philox_double2(time_step,
                       _mm256_add_epi32(
                           _mm256_add_epi32(
                               _mm256_set_epi32(block_offset_0, block_offset_0,
                                                block_offset_0, block_offset_0,
                                                block_offset_0, block_offset_0,
                                                block_offset_0, block_offset_0),
                               _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0)),
                           _mm256_set_epi32(ctr_0, ctr_0, ctr_0, ctr_0, ctr_0,
                                            ctr_0, ctr_0, ctr_0)),
                       block_offset_1 + ctr_1, block_offset_2 + ctr_2, 2, seed,
                       random_2_0, random_2_1);

        __m256d random_1_0;
        __m256d random_1_1;
        philox_double2(time_step,
                       _mm256_add_epi32(
                           _mm256_add_epi32(
                               _mm256_set_epi32(block_offset_0, block_offset_0,
                                                block_offset_0, block_offset_0,
                                                block_offset_0, block_offset_0,
                                                block_offset_0, block_offset_0),
                               _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0)),
                           _mm256_set_epi32(ctr_0, ctr_0, ctr_0, ctr_0, ctr_0,
                                            ctr_0, ctr_0, ctr_0)),
                       block_offset_1 + ctr_1, block_offset_2 + ctr_2, 1, seed,
                       random_1_0, random_1_1);

        __m256d random_0_0;
        __m256d random_0_1;
        philox_double2(time_step,
                       _mm256_add_epi32(
                           _mm256_add_epi32(
                               _mm256_set_epi32(block_offset_0, block_offset_0,
                                                block_offset_0, block_offset_0,
                                                block_offset_0, block_offset_0,
                                                block_offset_0, block_offset_0),
                               _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0)),
                           _mm256_set_epi32(ctr_0, ctr_0, ctr_0, ctr_0, ctr_0,
                                            ctr_0, ctr_0, ctr_0)),
                       block_offset_1 + ctr_1, block_offset_2 + ctr_2, 0, seed,
                       random_0_0, random_0_1);

        const __m256 xi_0 = _mm256_add_ps(xi_252, xi_266);
        const __m256 xi_1 = _mm256_add_ps(xi_0, xi_257);
        const __m256 xi_2 =
            _mm256_add_ps(_mm256_add_ps(xi_251, xi_262), xi_265);
        const __m256 xi_3 = _mm256_add_ps(xi_250, xi_253);
        const __m256 xi_4 = _mm256_add_ps(xi_263, xi_268);
        const __m256 xi_5 = _mm256_add_ps(xi_254, xi_267);
        const __m256 xi_6 = _mm256_add_ps(xi_260, xi_264);
        const __m256 xi_8 =
            _mm256_mul_ps(xi_268, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_9 = _mm256_add_ps(
            _mm256_mul_ps(xi_249, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)),
            xi_8);
        const __m256 xi_10 =
            _mm256_mul_ps(xi_260, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_11 =
            _mm256_mul_ps(xi_255, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_12 =
            _mm256_mul_ps(xi_263, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_13 = _mm256_add_ps(_mm256_add_ps(xi_10, xi_11), xi_12);
        const __m256 xi_14 =
            _mm256_mul_ps(xi_254, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_15 =
            _mm256_mul_ps(xi_256, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_16 = _mm256_add_ps(xi_14, xi_15);
        const __m256 xi_17 =
            _mm256_mul_ps(xi_267, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_18 =
            _mm256_mul_ps(xi_250, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_19 = _mm256_add_ps(xi_17, xi_18);
        const __m256 xi_20 =
            _mm256_mul_ps(xi_252, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_21 = _mm256_add_ps(xi_10, xi_20);
        const __m256 xi_22 =
            _mm256_mul_ps(xi_251, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_23 =
            _mm256_mul_ps(xi_264, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_24 = _mm256_add_ps(
            _mm256_add_ps(_mm256_add_ps(xi_17, xi_22), xi_23), xi_265);
        const __m256 xi_42 =
            _mm256_mul_ps(_mm256_set_ps(0.166666666666667, 0.166666666666667,
                                        0.166666666666667, 0.166666666666667,
                                        0.166666666666667, 0.166666666666667,
                                        0.166666666666667, 0.166666666666667),
                          _mm256_set_ps(xi_248, xi_248, xi_248, xi_248, xi_248,
                                        xi_248, xi_248, xi_248));
        const __m256 xi_50 =
            _mm256_mul_ps(_mm256_set_ps(0.166666666666667, 0.166666666666667,
                                        0.166666666666667, 0.166666666666667,
                                        0.166666666666667, 0.166666666666667,
                                        0.166666666666667, 0.166666666666667),
                          _mm256_set_ps(xi_269, xi_269, xi_269, xi_269, xi_269,
                                        xi_269, xi_269, xi_269));
        const __m256 xi_54 =
            _mm256_mul_ps(_mm256_set_ps(0.166666666666667, 0.166666666666667,
                                        0.166666666666667, 0.166666666666667,
                                        0.166666666666667, 0.166666666666667,
                                        0.166666666666667, 0.166666666666667),
                          _mm256_set_ps(xi_261, xi_261, xi_261, xi_261, xi_261,
                                        xi_261, xi_261, xi_261));
        const __m256 xi_57 =
            _mm256_mul_ps(_mm256_set_ps(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                          _mm256_set_ps(xi_248, xi_248, xi_248, xi_248, xi_248,
                                        xi_248, xi_248, xi_248));
        const __m256 xi_61 =
            _mm256_mul_ps(_mm256_set_ps(0.0833333333333333, 0.0833333333333333,
                                        0.0833333333333333, 0.0833333333333333,
                                        0.0833333333333333, 0.0833333333333333,
                                        0.0833333333333333, 0.0833333333333333),
                          _mm256_set_ps(xi_269, xi_269, xi_269, xi_269, xi_269,
                                        xi_269, xi_269, xi_269));
        const __m256 xi_65 =
            _mm256_mul_ps(_mm256_set_ps(0.0833333333333333, 0.0833333333333333,
                                        0.0833333333333333, 0.0833333333333333,
                                        0.0833333333333333, 0.0833333333333333,
                                        0.0833333333333333, 0.0833333333333333),
                          _mm256_set_ps(xi_248, xi_248, xi_248, xi_248, xi_248,
                                        xi_248, xi_248, xi_248));
        const __m256 xi_75 =
            _mm256_mul_ps(_mm256_set_ps(0.0833333333333333, 0.0833333333333333,
                                        0.0833333333333333, 0.0833333333333333,
                                        0.0833333333333333, 0.0833333333333333,
                                        0.0833333333333333, 0.0833333333333333),
                          _mm256_set_ps(xi_261, xi_261, xi_261, xi_261, xi_261,
                                        xi_261, xi_261, xi_261));
        const __m256 xi_93 =
            _mm256_mul_ps(xi_259, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_94 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_mul_ps(
                    _mm256_set_ps(3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0),
                    _mm256_set_ps(xi_253, xi_253, xi_253, xi_253, xi_253,
                                  xi_253, xi_253, xi_253)),
                _mm256_mul_ps(
                    _mm256_set_ps(3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0),
                    _mm256_set_ps(xi_264, xi_264, xi_264, xi_264, xi_264,
                                  xi_264, xi_264, xi_264))),
            _mm256_set_ps(xi_93, xi_93, xi_93, xi_93, xi_93, xi_93, xi_93,
                          xi_93));
        const __m256d xi_95 = _mm256_mul_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_mul_ps(
                                        _mm256_set_ps(-3.0, -3.0, -3.0, -3.0,
                                                      -3.0, -3.0, -3.0, -3.0),
                                        _mm256_set_ps(xi_250, xi_250, xi_250,
                                                      xi_250, xi_250, xi_250,
                                                      xi_250, xi_250)),
                                    _mm256_mul_ps(
                                        _mm256_set_ps(-3.0, -3.0, -3.0, -3.0,
                                                      -3.0, -3.0, -3.0, -3.0),
                                        _mm256_set_ps(xi_251, xi_251, xi_251,
                                                      xi_251, xi_251, xi_251,
                                                      xi_251, xi_251))),
                                _mm256_mul_ps(_mm256_set_ps(3.0, 3.0, 3.0, 3.0,
                                                            3.0, 3.0, 3.0, 3.0),
                                              _mm256_set_ps(xi_254, xi_254,
                                                            xi_254, xi_254,
                                                            xi_254, xi_254,
                                                            xi_254, xi_254))),
                            _mm256_mul_ps(_mm256_set_ps(3.0, 3.0, 3.0, 3.0, 3.0,
                                                        3.0, 3.0, 3.0),
                                          _mm256_set_ps(xi_262, xi_262, xi_262,
                                                        xi_262, xi_262, xi_262,
                                                        xi_262, xi_262))),
                        _mm256_mul_ps(_mm256_set_ps(-3.0, -3.0, -3.0, -3.0,
                                                    -3.0, -3.0, -3.0, -3.0),
                                      _mm256_set_ps(xi_265, xi_265, xi_265,
                                                    xi_265, xi_265, xi_265,
                                                    xi_265, xi_265))),
                    _mm256_mul_ps(_mm256_set_ps(-3.0, -3.0, -3.0, -3.0, -3.0,
                                                -3.0, -3.0, -3.0),
                                  _mm256_set_ps(xi_267, xi_267, xi_267, xi_267,
                                                xi_267, xi_267, xi_267,
                                                xi_267))),
                _mm256_set_ps(xi_94, xi_94, xi_94, xi_94, xi_94, xi_94, xi_94,
                              xi_94)),
            _mm256_set_ps(omega_even, omega_even, omega_even, omega_even,
                          omega_even, omega_even, omega_even, omega_even));
        const __m256 xi_96 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_mul_ps(
                        _mm256_set_ps(2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0),
                        _mm256_set_ps(xi_250, xi_250, xi_250, xi_250, xi_250,
                                      xi_250, xi_250, xi_250)),
                    _mm256_mul_ps(
                        _mm256_set_ps(2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0),
                        _mm256_set_ps(xi_251, xi_251, xi_251, xi_251, xi_251,
                                      xi_251, xi_251, xi_251))),
                _mm256_mul_ps(
                    _mm256_set_ps(2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0),
                    _mm256_set_ps(xi_265, xi_265, xi_265, xi_265, xi_265,
                                  xi_265, xi_265, xi_265))),
            _mm256_mul_ps(_mm256_set_ps(2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0),
                          _mm256_set_ps(xi_267, xi_267, xi_267, xi_267, xi_267,
                                        xi_267, xi_267, xi_267)));
        const __m256 xi_97 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_mul_ps(
                    _mm256_set_ps(5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0),
                    _mm256_set_ps(xi_257, xi_257, xi_257, xi_257, xi_257,
                                  xi_257, xi_257, xi_257)),
                _mm256_mul_ps(
                    _mm256_set_ps(5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0),
                    _mm256_set_ps(xi_263, xi_263, xi_263, xi_263, xi_263,
                                  xi_263, xi_263, xi_263))),
            _mm256_set_ps(xi_96, xi_96, xi_96, xi_96, xi_96, xi_96, xi_96,
                          xi_96));
        const __m256d xi_98 = _mm256_mul_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_add_ps(
                                        _mm256_mul_ps(
                                            _mm256_set_ps(-5.0, -5.0, -5.0,
                                                          -5.0, -5.0, -5.0,
                                                          -5.0, -5.0),
                                            _mm256_set_ps(xi_252, xi_252,
                                                          xi_252, xi_252,
                                                          xi_252, xi_252,
                                                          xi_252, xi_252)),
                                        _mm256_mul_ps(
                                            _mm256_set_ps(-2.0, -2.0, -2.0,
                                                          -2.0, -2.0, -2.0,
                                                          -2.0, -2.0),
                                            _mm256_set_ps(xi_254, xi_254,
                                                          xi_254, xi_254,
                                                          xi_254, xi_254,
                                                          xi_254, xi_254))),
                                    _mm256_mul_ps(
                                        _mm256_set_ps(-5.0, -5.0, -5.0, -5.0,
                                                      -5.0, -5.0, -5.0, -5.0),
                                        _mm256_set_ps(xi_255, xi_255, xi_255,
                                                      xi_255, xi_255, xi_255,
                                                      xi_255, xi_255))),
                                _mm256_mul_ps(
                                    _mm256_set_ps(-5.0, -5.0, -5.0, -5.0, -5.0,
                                                  -5.0, -5.0, -5.0),
                                    _mm256_set_ps(xi_260, xi_260, xi_260,
                                                  xi_260, xi_260, xi_260,
                                                  xi_260, xi_260))),
                            _mm256_mul_ps(_mm256_set_ps(-2.0, -2.0, -2.0, -2.0,
                                                        -2.0, -2.0, -2.0, -2.0),
                                          _mm256_set_ps(xi_262, xi_262, xi_262,
                                                        xi_262, xi_262, xi_262,
                                                        xi_262, xi_262))),
                        _mm256_mul_ps(_mm256_set_ps(-5.0, -5.0, -5.0, -5.0,
                                                    -5.0, -5.0, -5.0, -5.0),
                                      _mm256_set_ps(xi_266, xi_266, xi_266,
                                                    xi_266, xi_266, xi_266,
                                                    xi_266, xi_266))),
                    _mm256_set_ps(xi_94, xi_94, xi_94, xi_94, xi_94, xi_94,
                                  xi_94, xi_94)),
                _mm256_set_ps(xi_97, xi_97, xi_97, xi_97, xi_97, xi_97, xi_97,
                              xi_97)),
            _mm256_set_ps(omega_even, omega_even, omega_even, omega_even,
                          omega_even, omega_even, omega_even, omega_even));
        const __m256 xi_101 =
            _mm256_mul_ps(xi_265, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_102 = _mm256_add_ps(xi_101, xi_18);
        const __m256 xi_103 =
            _mm256_mul_ps(xi_258, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_106 =
            _mm256_mul_ps(xi_266, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_107 = _mm256_add_ps(
            _mm256_add_ps(_mm256_add_ps(xi_106, xi_11), xi_15), xi_21);
        const __m256 xi_109 =
            _mm256_mul_ps(_mm256_set_ps(2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0),
                          _mm256_set_ps(xi_255, xi_255, xi_255, xi_255, xi_255,
                                        xi_255, xi_255, xi_255));
        const __m256 xi_110 =
            _mm256_mul_ps(_mm256_set_ps(2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0),
                          _mm256_set_ps(xi_266, xi_266, xi_266, xi_266, xi_266,
                                        xi_266, xi_266, xi_266));
        const __m256 xi_111 = _mm256_add_ps(
            _mm256_mul_ps(_mm256_set_ps(2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0),
                          _mm256_set_ps(xi_252, xi_252, xi_252, xi_252, xi_252,
                                        xi_252, xi_252, xi_252)),
            _mm256_mul_ps(_mm256_set_ps(2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0),
                          _mm256_set_ps(xi_260, xi_260, xi_260, xi_260, xi_260,
                                        xi_260, xi_260, xi_260)));
        const __m256d xi_112 = _mm256_mul_ps(
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
                                                            _mm256_mul_ps(
                                                                _mm256_set_ps(
                                                                    -7.0, -7.0,
                                                                    -7.0, -7.0,
                                                                    -7.0, -7.0,
                                                                    -7.0, -7.0),
                                                                _mm256_set_ps(
                                                                    xi_249,
                                                                    xi_249,
                                                                    xi_249,
                                                                    xi_249,
                                                                    xi_249,
                                                                    xi_249,
                                                                    xi_249,
                                                                    xi_249)),
                                                            _mm256_mul_ps(
                                                                _mm256_set_ps(
                                                                    -4.0, -4.0,
                                                                    -4.0, -4.0,
                                                                    -4.0, -4.0,
                                                                    -4.0, -4.0),
                                                                _mm256_set_ps(
                                                                    xi_253,
                                                                    xi_253,
                                                                    xi_253,
                                                                    xi_253,
                                                                    xi_253,
                                                                    xi_253,
                                                                    xi_253,
                                                                    xi_253))),
                                                        _mm256_mul_ps(
                                                            _mm256_set_ps(
                                                                5.0, 5.0, 5.0,
                                                                5.0, 5.0, 5.0,
                                                                5.0, 5.0),
                                                            _mm256_set_ps(
                                                                xi_254, xi_254,
                                                                xi_254, xi_254,
                                                                xi_254, xi_254,
                                                                xi_254,
                                                                xi_254))),
                                                    _mm256_mul_ps(
                                                        _mm256_set_ps(
                                                            -7.0, -7.0, -7.0,
                                                            -7.0, -7.0, -7.0,
                                                            -7.0, -7.0),
                                                        _mm256_set_ps(
                                                            xi_256, xi_256,
                                                            xi_256, xi_256,
                                                            xi_256, xi_256,
                                                            xi_256, xi_256))),
                                                _mm256_mul_ps(
                                                    _mm256_set_ps(
                                                        -7.0, -7.0, -7.0, -7.0,
                                                        -7.0, -7.0, -7.0, -7.0),
                                                    _mm256_set_ps(
                                                        xi_258, xi_258, xi_258,
                                                        xi_258, xi_258, xi_258,
                                                        xi_258, xi_258))),
                                            _mm256_mul_ps(
                                                _mm256_set_ps(5.0, 5.0, 5.0,
                                                              5.0, 5.0, 5.0,
                                                              5.0, 5.0),
                                                _mm256_set_ps(xi_262, xi_262,
                                                              xi_262, xi_262,
                                                              xi_262, xi_262,
                                                              xi_262, xi_262))),
                                        _mm256_mul_ps(
                                            _mm256_set_ps(-4.0, -4.0, -4.0,
                                                          -4.0, -4.0, -4.0,
                                                          -4.0, -4.0),
                                            _mm256_set_ps(xi_264, xi_264,
                                                          xi_264, xi_264,
                                                          xi_264, xi_264,
                                                          xi_264, xi_264))),
                                    _mm256_mul_ps(
                                        _mm256_set_ps(-7.0, -7.0, -7.0, -7.0,
                                                      -7.0, -7.0, -7.0, -7.0),
                                        _mm256_set_ps(xi_268, xi_268, xi_268,
                                                      xi_268, xi_268, xi_268,
                                                      xi_268, xi_268))),
                                _mm256_set_ps(xi_109, xi_109, xi_109, xi_109,
                                              xi_109, xi_109, xi_109, xi_109)),
                            _mm256_set_ps(xi_110, xi_110, xi_110, xi_110,
                                          xi_110, xi_110, xi_110, xi_110)),
                        _mm256_set_ps(xi_111, xi_111, xi_111, xi_111, xi_111,
                                      xi_111, xi_111, xi_111)),
                    _mm256_set_ps(xi_93, xi_93, xi_93, xi_93, xi_93, xi_93,
                                  xi_93, xi_93)),
                _mm256_set_ps(xi_97, xi_97, xi_97, xi_97, xi_97, xi_97, xi_97,
                              xi_97)),
            _mm256_set_ps(omega_even, omega_even, omega_even, omega_even,
                          omega_even, omega_even, omega_even, omega_even));
        const __m256 xi_113 = _mm256_add_ps(xi_101, xi_250);
        const __m256 xi_114 = _mm256_add_ps(
            _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_113, xi_14), xi_22),
                          xi_262),
            xi_267);
        const __m256d xi_116 =
            _mm256_mul_ps(_mm256_set_ps(xi_114, xi_114, xi_114, xi_114, xi_114,
                                        xi_114, xi_114, xi_114),
                          _mm256_set_ps(xi_115, xi_115, xi_115, xi_115, xi_115,
                                        xi_115, xi_115, xi_115));
        const __m256 xi_118 = _mm256_add_ps(xi_103, xi_256);
        const __m256d xi_122 =
            _mm256_add_ps(random_5_1, _mm256_set_ps(-0.5, -0.5, -0.5, -0.5,
                                                    -0.5, -0.5, -0.5, -0.5));
        const __m256 xi_127 =
            _mm256_mul_ps(_mm256_set_ps(2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0),
                          _mm256_set_ps(xi_249, xi_249, xi_249, xi_249, xi_249,
                                        xi_249, xi_249, xi_249));
        const __m256 xi_128 =
            _mm256_mul_ps(_mm256_set_ps(2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0),
                          _mm256_set_ps(xi_256, xi_256, xi_256, xi_256, xi_256,
                                        xi_256, xi_256, xi_256));
        const __m256 xi_129 = _mm256_add_ps(
            _mm256_mul_ps(
                _mm256_set_ps(-2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0),
                _mm256_set_ps(xi_258, xi_258, xi_258, xi_258, xi_258, xi_258,
                              xi_258, xi_258)),
            _mm256_mul_ps(_mm256_set_ps(2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0),
                          _mm256_set_ps(xi_268, xi_268, xi_268, xi_268, xi_268,
                                        xi_268, xi_268, xi_268)));
        const __m256 xi_130 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_mul_ps(
                                xi_127, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                            xi_128),
                        xi_129),
                    xi_14),
                xi_19),
            xi_2);
        const __m256d xi_132 =
            _mm256_mul_ps(_mm256_set_ps(xi_130, xi_130, xi_130, xi_130, xi_130,
                                        xi_130, xi_130, xi_130),
                          _mm256_set_ps(xi_131, xi_131, xi_131, xi_131, xi_131,
                                        xi_131, xi_131, xi_131));
        const __m256d xi_133 =
            _mm256_add_ps(random_3_0, _mm256_set_ps(-0.5, -0.5, -0.5, -0.5,
                                                    -0.5, -0.5, -0.5, -0.5));
        const __m256d xi_138 =
            _mm256_add_ps(random_0_1, _mm256_set_ps(-0.5, -0.5, -0.5, -0.5,
                                                    -0.5, -0.5, -0.5, -0.5));
        const __m256 xi_142 = _mm256_add_ps(xi_255, xi_260);
        const __m256 xi_156 = _mm256_add_ps(xi_106, xi_255);
        const __m256 xi_157 = _mm256_add_ps(
            _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_12, xi_156), xi_20),
                          xi_257),
            xi_260);
        const __m256d xi_158 =
            _mm256_mul_ps(_mm256_set_ps(xi_115, xi_115, xi_115, xi_115, xi_115,
                                        xi_115, xi_115, xi_115),
                          _mm256_set_ps(xi_157, xi_157, xi_157, xi_157, xi_157,
                                        xi_157, xi_157, xi_157));
        const __m256d xi_159 =
            _mm256_add_ps(random_4_1, _mm256_set_ps(-0.5, -0.5, -0.5, -0.5,
                                                    -0.5, -0.5, -0.5, -0.5));
        const __m256 xi_161 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_mul_ps(xi_128,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        xi_1),
                    xi_127),
                xi_129),
            xi_13);
        const __m256d xi_162 =
            _mm256_mul_ps(_mm256_set_ps(xi_131, xi_131, xi_131, xi_131, xi_131,
                                        xi_131, xi_131, xi_131),
                          _mm256_set_ps(xi_161, xi_161, xi_161, xi_161, xi_161,
                                        xi_161, xi_161, xi_161));
        const __m256d xi_163 =
            _mm256_add_ps(random_4_0, _mm256_set_ps(-0.5, -0.5, -0.5, -0.5,
                                                    -0.5, -0.5, -0.5, -0.5));
        const __m256 xi_168 = _mm256_add_ps(xi_251, xi_267);
        const __m256 xi_169 = _mm256_add_ps(
            _mm256_add_ps(_mm256_add_ps(xi_102, xi_168), xi_23), xi_253);
        const __m256d xi_170 =
            _mm256_mul_ps(_mm256_set_ps(xi_115, xi_115, xi_115, xi_115, xi_115,
                                        xi_115, xi_115, xi_115),
                          _mm256_set_ps(xi_169, xi_169, xi_169, xi_169, xi_169,
                                        xi_169, xi_169, xi_169));
        const __m256d xi_173 =
            _mm256_add_ps(random_5_0, _mm256_set_ps(-0.5, -0.5, -0.5, -0.5,
                                                    -0.5, -0.5, -0.5, -0.5));
        const __m256 xi_175 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_mul_ps(xi_109,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_110,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0))),
                    xi_111),
                xi_24),
            xi_3);
        const __m256d xi_176 =
            _mm256_mul_ps(_mm256_set_ps(xi_131, xi_131, xi_131, xi_131, xi_131,
                                        xi_131, xi_131, xi_131),
                          _mm256_set_ps(xi_175, xi_175, xi_175, xi_175, xi_175,
                                        xi_175, xi_175, xi_175));
        const __m256d xi_177 =
            _mm256_add_ps(random_3_1, _mm256_set_ps(-0.5, -0.5, -0.5, -0.5,
                                                    -0.5, -0.5, -0.5, -0.5));
        const __m256d xi_184 = _mm256_mul_ps(
            xi_112, _mm256_set_ps(0.0138888888888889, 0.0138888888888889,
                                  0.0138888888888889, 0.0138888888888889,
                                  0.0138888888888889, 0.0138888888888889,
                                  0.0138888888888889, 0.0138888888888889));
        const __m256d xi_205 = _mm256_mul_ps(
            xi_98, _mm256_set_ps(-0.00714285714285714, -0.00714285714285714,
                                 -0.00714285714285714, -0.00714285714285714,
                                 -0.00714285714285714, -0.00714285714285714,
                                 -0.00714285714285714, -0.00714285714285714));
        const __m256d xi_207 =
            _mm256_mul_ps(xi_95, _mm256_set_ps(0.025, 0.025, 0.025, 0.025,
                                               0.025, 0.025, 0.025, 0.025));
        const __m256d xi_212 =
            _mm256_mul_ps(_mm256_set_ps(xi_175, xi_175, xi_175, xi_175, xi_175,
                                        xi_175, xi_175, xi_175),
                          _mm256_set_ps(xi_211, xi_211, xi_211, xi_211, xi_211,
                                        xi_211, xi_211, xi_211));
        const __m256d xi_214 =
            _mm256_mul_ps(_mm256_set_ps(xi_169, xi_169, xi_169, xi_169, xi_169,
                                        xi_169, xi_169, xi_169),
                          _mm256_set_ps(xi_213, xi_213, xi_213, xi_213, xi_213,
                                        xi_213, xi_213, xi_213));
        const __m256d xi_223 =
            _mm256_mul_ps(_mm256_set_ps(xi_130, xi_130, xi_130, xi_130, xi_130,
                                        xi_130, xi_130, xi_130),
                          _mm256_set_ps(xi_211, xi_211, xi_211, xi_211, xi_211,
                                        xi_211, xi_211, xi_211));
        const __m256d xi_224 =
            _mm256_mul_ps(_mm256_set_ps(xi_114, xi_114, xi_114, xi_114, xi_114,
                                        xi_114, xi_114, xi_114),
                          _mm256_set_ps(xi_213, xi_213, xi_213, xi_213, xi_213,
                                        xi_213, xi_213, xi_213));
        const __m256d xi_232 = _mm256_mul_ps(
            xi_98, _mm256_set_ps(0.0178571428571429, 0.0178571428571429,
                                 0.0178571428571429, 0.0178571428571429,
                                 0.0178571428571429, 0.0178571428571429,
                                 0.0178571428571429, 0.0178571428571429));
        const __m256d xi_238 =
            _mm256_mul_ps(_mm256_set_ps(xi_157, xi_157, xi_157, xi_157, xi_157,
                                        xi_157, xi_157, xi_157),
                          _mm256_set_ps(xi_213, xi_213, xi_213, xi_213, xi_213,
                                        xi_213, xi_213, xi_213));
        const __m256d xi_239 =
            _mm256_mul_ps(_mm256_set_ps(xi_161, xi_161, xi_161, xi_161, xi_161,
                                        xi_161, xi_161, xi_161),
                          _mm256_set_ps(xi_211, xi_211, xi_211, xi_211, xi_211,
                                        xi_211, xi_211, xi_211));
        const __m256 vel0Term =
            _mm256_add_ps(_mm256_add_ps(xi_1, xi_256), xi_258);
        const __m256 vel1Term = _mm256_add_ps(xi_2, xi_249);
        const __m256 vel2Term = _mm256_add_ps(xi_255, xi_3);
        const __m256 rho = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_add_ps(vel0Term, vel1Term),
                                      vel2Term),
                        xi_259),
                    xi_4),
                xi_5),
            xi_6);
        const __m256 xi_7 = _mm256_div_ps(
            _mm256_set_ps(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0), rho);
        const __m256d xi_86 = _mm256_mul_ps(
            _mm256_set_ps(kT, kT, kT, kT, kT, kT, kT, kT),
            _mm256_set_ps(rho, rho, rho, rho, rho, rho, rho, rho));
        const __m256d xi_87 = _mm256_sqrt_ps(_mm256_mul_ps(
            xi_86,
            _mm256_set_ps(-((-omega_even + 1.0) * (-omega_even + 1.0)) + 1.0,
                          -((-omega_even + 1.0) * (-omega_even + 1.0)) + 1.0,
                          -((-omega_even + 1.0) * (-omega_even + 1.0)) + 1.0,
                          -((-omega_even + 1.0) * (-omega_even + 1.0)) + 1.0,
                          -((-omega_even + 1.0) * (-omega_even + 1.0)) + 1.0,
                          -((-omega_even + 1.0) * (-omega_even + 1.0)) + 1.0,
                          -((-omega_even + 1.0) * (-omega_even + 1.0)) + 1.0,
                          -((-omega_even + 1.0) * (-omega_even + 1.0)) + 1.0)));
        const __m256d xi_88 = _mm256_mul_ps(
            _mm256_mul_ps(xi_87,
                          _mm256_add_ps(random_6_0,
                                        _mm256_set_ps(-0.5, -0.5, -0.5, -0.5,
                                                      -0.5, -0.5, -0.5, -0.5))),
            _mm256_set_ps(3.7416573867739413, 3.7416573867739413,
                          3.7416573867739413, 3.7416573867739413,
                          3.7416573867739413, 3.7416573867739413,
                          3.7416573867739413, 3.7416573867739413));
        const __m256d xi_89 = _mm256_mul_ps(
            _mm256_mul_ps(xi_87,
                          _mm256_add_ps(random_7_0,
                                        _mm256_set_ps(-0.5, -0.5, -0.5, -0.5,
                                                      -0.5, -0.5, -0.5, -0.5))),
            _mm256_set_ps(5.4772255750516612, 5.4772255750516612,
                          5.4772255750516612, 5.4772255750516612,
                          5.4772255750516612, 5.4772255750516612,
                          5.4772255750516612, 5.4772255750516612));
        const __m256d xi_91 = _mm256_mul_ps(
            _mm256_mul_ps(
                _mm256_sqrt_ps(_mm256_mul_ps(
                    xi_86,
                    _mm256_set_ps(-((xi_25 + 1.0) * (xi_25 + 1.0)) + 1.0,
                                  -((xi_25 + 1.0) * (xi_25 + 1.0)) + 1.0,
                                  -((xi_25 + 1.0) * (xi_25 + 1.0)) + 1.0,
                                  -((xi_25 + 1.0) * (xi_25 + 1.0)) + 1.0,
                                  -((xi_25 + 1.0) * (xi_25 + 1.0)) + 1.0,
                                  -((xi_25 + 1.0) * (xi_25 + 1.0)) + 1.0,
                                  -((xi_25 + 1.0) * (xi_25 + 1.0)) + 1.0,
                                  -((xi_25 + 1.0) * (xi_25 + 1.0)) + 1.0))),
                _mm256_add_ps(random_2_1,
                              _mm256_set_ps(-0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
                                            -0.5, -0.5))),
            _mm256_set_ps(xi_90, xi_90, xi_90, xi_90, xi_90, xi_90, xi_90,
                          xi_90));
        const __m256d xi_92 = _mm256_mul_ps(
            _mm256_mul_ps(xi_87,
                          _mm256_add_ps(random_6_1,
                                        _mm256_set_ps(-0.5, -0.5, -0.5, -0.5,
                                                      -0.5, -0.5, -0.5, -0.5))),
            _mm256_set_ps(8.3666002653407556, 8.3666002653407556,
                          8.3666002653407556, 8.3666002653407556,
                          8.3666002653407556, 8.3666002653407556,
                          8.3666002653407556, 8.3666002653407556));
        const __m256d xi_123 = _mm256_sqrt_ps(_mm256_mul_ps(
            xi_86,
            _mm256_set_ps(-((-omega_odd + 1.0) * (-omega_odd + 1.0)) + 1.0,
                          -((-omega_odd + 1.0) * (-omega_odd + 1.0)) + 1.0,
                          -((-omega_odd + 1.0) * (-omega_odd + 1.0)) + 1.0,
                          -((-omega_odd + 1.0) * (-omega_odd + 1.0)) + 1.0,
                          -((-omega_odd + 1.0) * (-omega_odd + 1.0)) + 1.0,
                          -((-omega_odd + 1.0) * (-omega_odd + 1.0)) + 1.0,
                          -((-omega_odd + 1.0) * (-omega_odd + 1.0)) + 1.0,
                          -((-omega_odd + 1.0) * (-omega_odd + 1.0)) + 1.0)));
        const __m256d xi_124 = _mm256_mul_ps(
            xi_123, _mm256_set_ps(1.4142135623730951, 1.4142135623730951,
                                  1.4142135623730951, 1.4142135623730951,
                                  1.4142135623730951, 1.4142135623730951,
                                  1.4142135623730951, 1.4142135623730951));
        const __m256d xi_125 = _mm256_mul_ps(
            xi_124, _mm256_set_ps(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5));
        const __m256d xi_126 = _mm256_mul_ps(xi_122, xi_125);
        const __m256d xi_134 =
            _mm256_mul_ps(xi_123, _mm256_set_ps(xi_90, xi_90, xi_90, xi_90,
                                                xi_90, xi_90, xi_90, xi_90));
        const __m256d xi_135 = _mm256_mul_ps(
            xi_134, _mm256_set_ps(0.166666666666667, 0.166666666666667,
                                  0.166666666666667, 0.166666666666667,
                                  0.166666666666667, 0.166666666666667,
                                  0.166666666666667, 0.166666666666667));
        const __m256d xi_136 = _mm256_mul_ps(xi_133, xi_135);
        const __m256d xi_137 = _mm256_add_ps(
            _mm256_mul_ps(xi_132, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)),
            _mm256_mul_ps(xi_136, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)));
        const __m256d xi_139 = _mm256_sqrt_ps(_mm256_mul_ps(
            xi_86, _mm256_set_ps(-((xi_36 + 1.0) * (xi_36 + 1.0)) + 1.0,
                                 -((xi_36 + 1.0) * (xi_36 + 1.0)) + 1.0,
                                 -((xi_36 + 1.0) * (xi_36 + 1.0)) + 1.0,
                                 -((xi_36 + 1.0) * (xi_36 + 1.0)) + 1.0,
                                 -((xi_36 + 1.0) * (xi_36 + 1.0)) + 1.0,
                                 -((xi_36 + 1.0) * (xi_36 + 1.0)) + 1.0,
                                 -((xi_36 + 1.0) * (xi_36 + 1.0)) + 1.0,
                                 -((xi_36 + 1.0) * (xi_36 + 1.0)) + 1.0)));
        const __m256d xi_140 = _mm256_mul_ps(
            xi_139, _mm256_set_ps(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5));
        const __m256d xi_141 = _mm256_mul_ps(xi_138, xi_140);
        const __m256d xi_146 = _mm256_add_ps(
            _mm256_mul_ps(
                xi_112,
                _mm256_set_ps(-0.0198412698412698, -0.0198412698412698,
                              -0.0198412698412698, -0.0198412698412698,
                              -0.0198412698412698, -0.0198412698412698,
                              -0.0198412698412698, -0.0198412698412698)),
            _mm256_mul_ps(
                xi_88, _mm256_set_ps(-0.119047619047619, -0.119047619047619,
                                     -0.119047619047619, -0.119047619047619,
                                     -0.119047619047619, -0.119047619047619,
                                     -0.119047619047619, -0.119047619047619)));
        const __m256d xi_148 = _mm256_mul_ps(
            _mm256_mul_ps(xi_139,
                          _mm256_add_ps(random_0_0,
                                        _mm256_set_ps(-0.5, -0.5, -0.5, -0.5,
                                                      -0.5, -0.5, -0.5, -0.5))),
            _mm256_set_ps(1.7320508075688772, 1.7320508075688772,
                          1.7320508075688772, 1.7320508075688772,
                          1.7320508075688772, 1.7320508075688772,
                          1.7320508075688772, 1.7320508075688772));
        const __m256d xi_152 = _mm256_add_ps(xi_132, xi_136);
        const __m256d xi_160 = _mm256_mul_ps(xi_125, xi_159);
        const __m256d xi_164 = _mm256_mul_ps(xi_135, xi_163);
        const __m256d xi_165 = _mm256_add_ps(xi_162, xi_164);
        const __m256d xi_167 = _mm256_add_ps(
            _mm256_mul_ps(xi_162, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)),
            _mm256_mul_ps(xi_164, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)));
        const __m256d xi_174 = _mm256_mul_ps(xi_125, xi_173);
        const __m256d xi_178 = _mm256_mul_ps(xi_135, xi_177);
        const __m256d xi_179 = _mm256_add_ps(
            _mm256_mul_ps(xi_176, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)),
            _mm256_mul_ps(xi_178, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)));
        const __m256d xi_181 = _mm256_add_ps(xi_176, xi_178);
        const __m256d xi_182 = _mm256_mul_ps(
            _mm256_mul_ps(xi_138, xi_139),
            _mm256_set_ps(0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25));
        const __m256d xi_185 = _mm256_mul_ps(
            xi_88, _mm256_set_ps(0.0833333333333333, 0.0833333333333333,
                                 0.0833333333333333, 0.0833333333333333,
                                 0.0833333333333333, 0.0833333333333333,
                                 0.0833333333333333, 0.0833333333333333));
        const __m256d xi_195 = _mm256_mul_ps(
            xi_140,
            _mm256_add_ps(random_1_0, _mm256_set_ps(-0.5, -0.5, -0.5, -0.5,
                                                    -0.5, -0.5, -0.5, -0.5)));
        const __m256d xi_204 = _mm256_mul_ps(
            xi_140,
            _mm256_add_ps(random_2_0, _mm256_set_ps(-0.5, -0.5, -0.5, -0.5,
                                                    -0.5, -0.5, -0.5, -0.5)));
        const __m256d xi_208 = _mm256_mul_ps(
            xi_92, _mm256_set_ps(-0.0142857142857143, -0.0142857142857143,
                                 -0.0142857142857143, -0.0142857142857143,
                                 -0.0142857142857143, -0.0142857142857143,
                                 -0.0142857142857143, -0.0142857142857143));
        const __m256d xi_209 =
            _mm256_mul_ps(xi_89, _mm256_set_ps(0.05, 0.05, 0.05, 0.05, 0.05,
                                               0.05, 0.05, 0.05));
        const __m256d xi_215 = _mm256_mul_ps(
            xi_134, _mm256_set_ps(0.0833333333333333, 0.0833333333333333,
                                  0.0833333333333333, 0.0833333333333333,
                                  0.0833333333333333, 0.0833333333333333,
                                  0.0833333333333333, 0.0833333333333333));
        const __m256d xi_216 = _mm256_mul_ps(xi_177, xi_215);
        const __m256d xi_217 =
            _mm256_mul_ps(xi_124, _mm256_set_ps(0.25, 0.25, 0.25, 0.25, 0.25,
                                                0.25, 0.25, 0.25));
        const __m256d xi_218 = _mm256_mul_ps(xi_173, xi_217);
        const __m256d xi_220 = _mm256_add_ps(
            _mm256_mul_ps(
                xi_112,
                _mm256_set_ps(-0.00396825396825397, -0.00396825396825397,
                              -0.00396825396825397, -0.00396825396825397,
                              -0.00396825396825397, -0.00396825396825397,
                              -0.00396825396825397, -0.00396825396825397)),
            _mm256_mul_ps(
                xi_88,
                _mm256_set_ps(-0.0238095238095238, -0.0238095238095238,
                              -0.0238095238095238, -0.0238095238095238,
                              -0.0238095238095238, -0.0238095238095238,
                              -0.0238095238095238, -0.0238095238095238)));
        const __m256d xi_225 = _mm256_mul_ps(xi_133, xi_215);
        const __m256d xi_226 = _mm256_mul_ps(xi_122, xi_217);
        const __m256d xi_230 =
            _mm256_mul_ps(xi_182, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256d xi_233 = _mm256_mul_ps(
            xi_92, _mm256_set_ps(0.0357142857142857, 0.0357142857142857,
                                 0.0357142857142857, 0.0357142857142857,
                                 0.0357142857142857, 0.0357142857142857,
                                 0.0357142857142857, 0.0357142857142857));
        const __m256d xi_235 = _mm256_mul_ps(
            xi_140,
            _mm256_add_ps(random_1_1, _mm256_set_ps(-0.5, -0.5, -0.5, -0.5,
                                                    -0.5, -0.5, -0.5, -0.5)));
        const __m256d xi_240 = _mm256_mul_ps(xi_159, xi_217);
        const __m256d xi_241 = _mm256_mul_ps(xi_163, xi_215);
        const __m256 u_0 = _mm256_mul_ps(
            xi_7, _mm256_add_ps(_mm256_add_ps(vel0Term, xi_13), xi_9));
        const __m256 xi_26 = _mm256_mul_ps(u_0, xi_269);
        const __m256 xi_27 =
            _mm256_mul_ps(_mm256_set_ps(0.333333333333333, 0.333333333333333,
                                        0.333333333333333, 0.333333333333333,
                                        0.333333333333333, 0.333333333333333,
                                        0.333333333333333, 0.333333333333333),
                          _mm256_set_ps(xi_26, xi_26, xi_26, xi_26, xi_26,
                                        xi_26, xi_26, xi_26));
        const __m256 xi_33 =
            _mm256_mul_ps(xi_27, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                               -1.0, -1.0, -1.0));
        const __m256 xi_99 = _mm256_mul_ps(rho, (_mm256_mul_ps(u_0, u_0)));
        const __m256 xi_153 = _mm256_mul_ps(rho, u_0);
        const __m256 xi_154 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_mul_ps(vel0Term,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        xi_142),
                    xi_153),
                xi_249),
            xi_4);
        const __m256 xi_155 = _mm256_mul_ps(
            xi_154, _mm256_set_ps(xi_120, xi_120, xi_120, xi_120, xi_120,
                                  xi_120, xi_120, xi_120));
        const __m256 xi_191 = _mm256_mul_ps(
            xi_154, _mm256_set_ps(xi_186, xi_186, xi_186, xi_186, xi_186,
                                  xi_186, xi_186, xi_186));
        const __m256 u_1 = _mm256_mul_ps(
            xi_7, _mm256_add_ps(
                      _mm256_add_ps(
                          _mm256_add_ps(_mm256_add_ps(vel1Term, xi_16), xi_19),
                          xi_258),
                      xi_8));
        const __m256 xi_28 = _mm256_mul_ps(u_1, xi_248);
        const __m256 xi_29 =
            _mm256_mul_ps(_mm256_set_ps(0.333333333333333, 0.333333333333333,
                                        0.333333333333333, 0.333333333333333,
                                        0.333333333333333, 0.333333333333333,
                                        0.333333333333333, 0.333333333333333),
                          _mm256_set_ps(xi_28, xi_28, xi_28, xi_28, xi_28,
                                        xi_28, xi_28, xi_28));
        const __m256 xi_34 =
            _mm256_mul_ps(xi_29, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                               -1.0, -1.0, -1.0));
        const __m256 xi_56 = _mm256_mul_ps(
            _mm256_set_ps(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
            _mm256_set_ps(u_1, u_1, u_1, u_1, u_1, u_1, u_1, u_1));
        const __m256d xi_59 = _mm256_mul_ps(
            _mm256_set_ps(_mm256_add_ps(_mm256_mul_ps(u_0, xi_57),
                                        _mm256_mul_ps(xi_269, xi_56)),
                          _mm256_add_ps(_mm256_mul_ps(u_0, xi_57),
                                        _mm256_mul_ps(xi_269, xi_56)),
                          _mm256_add_ps(_mm256_mul_ps(u_0, xi_57),
                                        _mm256_mul_ps(xi_269, xi_56)),
                          _mm256_add_ps(_mm256_mul_ps(u_0, xi_57),
                                        _mm256_mul_ps(xi_269, xi_56)),
                          _mm256_add_ps(_mm256_mul_ps(u_0, xi_57),
                                        _mm256_mul_ps(xi_269, xi_56)),
                          _mm256_add_ps(_mm256_mul_ps(u_0, xi_57),
                                        _mm256_mul_ps(xi_269, xi_56)),
                          _mm256_add_ps(_mm256_mul_ps(u_0, xi_57),
                                        _mm256_mul_ps(xi_269, xi_56)),
                          _mm256_add_ps(_mm256_mul_ps(u_0, xi_57),
                                        _mm256_mul_ps(xi_269, xi_56))),
            _mm256_set_ps(xi_58, xi_58, xi_58, xi_58, xi_58, xi_58, xi_58,
                          xi_58));
        const __m256d xi_60 =
            _mm256_mul_ps(xi_59, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                               -1.0, -1.0, -1.0));
        const __m256 xi_104 = _mm256_mul_ps(rho, (_mm256_mul_ps(u_1, u_1)));
        const __m256 xi_105 =
            _mm256_add_ps(_mm256_add_ps(xi_103, xi_104), xi_9);
        const __m256 xi_117 = _mm256_mul_ps(rho, u_1);
        const __m256 xi_119 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_mul_ps(vel1Term,
                                          _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                        -1.0, -1.0, -1.0,
                                                        -1.0)),
                            xi_117),
                        xi_118),
                    xi_250),
                xi_268),
            xi_5);
        const __m256 xi_121 = _mm256_mul_ps(
            xi_119, _mm256_set_ps(xi_120, xi_120, xi_120, xi_120, xi_120,
                                  xi_120, xi_120, xi_120));
        const __m256 xi_187 = _mm256_mul_ps(
            xi_119, _mm256_set_ps(xi_186, xi_186, xi_186, xi_186, xi_186,
                                  xi_186, xi_186, xi_186));
        const __m256d xi_197 = _mm256_mul_ps(
            _mm256_set_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_mul_ps(u_0, xi_117), xi_118),
                        xi_249),
                    xi_8),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_mul_ps(u_0, xi_117), xi_118),
                        xi_249),
                    xi_8),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_mul_ps(u_0, xi_117), xi_118),
                        xi_249),
                    xi_8),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_mul_ps(u_0, xi_117), xi_118),
                        xi_249),
                    xi_8),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_mul_ps(u_0, xi_117), xi_118),
                        xi_249),
                    xi_8),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_mul_ps(u_0, xi_117), xi_118),
                        xi_249),
                    xi_8),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_mul_ps(u_0, xi_117), xi_118),
                        xi_249),
                    xi_8),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_mul_ps(u_0, xi_117), xi_118),
                        xi_249),
                    xi_8)),
            _mm256_set_ps(xi_196, xi_196, xi_196, xi_196, xi_196, xi_196,
                          xi_196, xi_196));
        const __m256d xi_198 = _mm256_add_ps(
            _mm256_mul_ps(xi_195, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)),
            _mm256_mul_ps(xi_197, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)));
        const __m256d xi_199 = _mm256_add_ps(xi_195, xi_197);
        const __m256 u_2 = _mm256_mul_ps(
            xi_7,
            _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(vel2Term, xi_21), xi_24),
                          xi_266));
        const __m256 xi_30 = _mm256_mul_ps(u_2, xi_261);
        const __m256 xi_31 =
            _mm256_mul_ps(_mm256_set_ps(0.333333333333333, 0.333333333333333,
                                        0.333333333333333, 0.333333333333333,
                                        0.333333333333333, 0.333333333333333,
                                        0.333333333333333, 0.333333333333333),
                          _mm256_set_ps(xi_30, xi_30, xi_30, xi_30, xi_30,
                                        xi_30, xi_30, xi_30));
        const __m256d xi_32 = _mm256_mul_ps(
            _mm256_set_ps(xi_25 + 2.0, xi_25 + 2.0, xi_25 + 2.0, xi_25 + 2.0,
                          xi_25 + 2.0, xi_25 + 2.0, xi_25 + 2.0, xi_25 + 2.0),
            _mm256_set_ps(_mm256_add_ps(_mm256_add_ps(xi_27, xi_29), xi_31),
                          _mm256_add_ps(_mm256_add_ps(xi_27, xi_29), xi_31),
                          _mm256_add_ps(_mm256_add_ps(xi_27, xi_29), xi_31),
                          _mm256_add_ps(_mm256_add_ps(xi_27, xi_29), xi_31),
                          _mm256_add_ps(_mm256_add_ps(xi_27, xi_29), xi_31),
                          _mm256_add_ps(_mm256_add_ps(xi_27, xi_29), xi_31),
                          _mm256_add_ps(_mm256_add_ps(xi_27, xi_29), xi_31),
                          _mm256_add_ps(_mm256_add_ps(xi_27, xi_29), xi_31)));
        const __m256 xi_35 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_mul_ps(
                    _mm256_set_ps(0.666666666666667, 0.666666666666667,
                                  0.666666666666667, 0.666666666666667,
                                  0.666666666666667, 0.666666666666667,
                                  0.666666666666667, 0.666666666666667),
                    _mm256_set_ps(xi_30, xi_30, xi_30, xi_30, xi_30, xi_30,
                                  xi_30, xi_30)),
                _mm256_set_ps(xi_33, xi_33, xi_33, xi_33, xi_33, xi_33, xi_33,
                              xi_33)),
            _mm256_set_ps(xi_34, xi_34, xi_34, xi_34, xi_34, xi_34, xi_34,
                          xi_34));
        const __m256 xi_39 =
            _mm256_mul_ps(xi_31, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                               -1.0, -1.0, -1.0));
        const __m256 xi_40 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_mul_ps(
                    _mm256_set_ps(0.666666666666667, 0.666666666666667,
                                  0.666666666666667, 0.666666666666667,
                                  0.666666666666667, 0.666666666666667,
                                  0.666666666666667, 0.666666666666667),
                    _mm256_set_ps(xi_28, xi_28, xi_28, xi_28, xi_28, xi_28,
                                  xi_28, xi_28)),
                _mm256_set_ps(xi_33, xi_33, xi_33, xi_33, xi_33, xi_33, xi_33,
                              xi_33)),
            _mm256_set_ps(xi_39, xi_39, xi_39, xi_39, xi_39, xi_39, xi_39,
                          xi_39));
        const __m256 xi_41 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_mul_ps(
                    _mm256_set_ps(0.666666666666667, 0.666666666666667,
                                  0.666666666666667, 0.666666666666667,
                                  0.666666666666667, 0.666666666666667,
                                  0.666666666666667, 0.666666666666667),
                    _mm256_set_ps(xi_26, xi_26, xi_26, xi_26, xi_26, xi_26,
                                  xi_26, xi_26)),
                _mm256_set_ps(xi_34, xi_34, xi_34, xi_34, xi_34, xi_34, xi_34,
                              xi_34)),
            _mm256_set_ps(xi_39, xi_39, xi_39, xi_39, xi_39, xi_39, xi_39,
                          xi_39));
        const __m256d xi_44 =
            _mm256_mul_ps(_mm256_set_ps(xi_35, xi_35, xi_35, xi_35, xi_35,
                                        xi_35, xi_35, xi_35),
                          _mm256_set_ps(xi_43, xi_43, xi_43, xi_43, xi_43,
                                        xi_43, xi_43, xi_43));
        const __m256d xi_45 =
            _mm256_mul_ps(xi_44, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                               -1.0, -1.0, -1.0));
        const __m256d xi_46 =
            _mm256_mul_ps(_mm256_set_ps(xi_41, xi_41, xi_41, xi_41, xi_41,
                                        xi_41, xi_41, xi_41),
                          _mm256_set_ps(xi_43, xi_43, xi_43, xi_43, xi_43,
                                        xi_43, xi_43, xi_43));
        const __m256d xi_47 =
            _mm256_mul_ps(xi_46, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                               -1.0, -1.0, -1.0));
        const __m256d xi_49 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_mul_ps(_mm256_set_ps(xi_40, xi_40, xi_40, xi_40, xi_40,
                                            xi_40, xi_40, xi_40),
                              _mm256_set_ps(xi_48, xi_48, xi_48, xi_48, xi_48,
                                            xi_48, xi_48, xi_48)),
                xi_45),
            xi_47);
        const __m256d xi_51 =
            _mm256_mul_ps(_mm256_set_ps(xi_40, xi_40, xi_40, xi_40, xi_40,
                                        xi_40, xi_40, xi_40),
                          _mm256_set_ps(xi_43, xi_43, xi_43, xi_43, xi_43,
                                        xi_43, xi_43, xi_43));
        const __m256d xi_52 =
            _mm256_mul_ps(xi_51, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                               -1.0, -1.0, -1.0));
        const __m256d xi_53 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_mul_ps(_mm256_set_ps(xi_41, xi_41, xi_41, xi_41, xi_41,
                                            xi_41, xi_41, xi_41),
                              _mm256_set_ps(xi_48, xi_48, xi_48, xi_48, xi_48,
                                            xi_48, xi_48, xi_48)),
                xi_45),
            xi_52);
        const __m256d xi_55 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_mul_ps(_mm256_set_ps(xi_35, xi_35, xi_35, xi_35, xi_35,
                                            xi_35, xi_35, xi_35),
                              _mm256_set_ps(xi_48, xi_48, xi_48, xi_48, xi_48,
                                            xi_48, xi_48, xi_48)),
                xi_47),
            xi_52);
        const __m256d xi_62 = _mm256_add_ps(
            xi_46,
            _mm256_set_ps(
                _mm256_mul_ps(xi_61, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_61, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_61, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_61, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_61, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_61, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_61, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_61, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0))));
        const __m256d xi_64 = _mm256_mul_ps(
            _mm256_mul_ps(
                _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0),
                _mm256_set_ps(xi_35, xi_35, xi_35, xi_35, xi_35, xi_35, xi_35,
                              xi_35)),
            _mm256_set_ps(xi_63, xi_63, xi_63, xi_63, xi_63, xi_63, xi_63,
                          xi_63));
        const __m256d xi_66 =
            _mm256_mul_ps(xi_32, _mm256_set_ps(0.125, 0.125, 0.125, 0.125,
                                               0.125, 0.125, 0.125, 0.125));
        const __m256d xi_67 = _mm256_add_ps(xi_51, xi_66);
        const __m256d xi_68 =
            _mm256_add_ps(xi_67, _mm256_set_ps(xi_65, xi_65, xi_65, xi_65,
                                               xi_65, xi_65, xi_65, xi_65));
        const __m256d xi_69 = _mm256_add_ps(xi_64, xi_68);
        const __m256d xi_70 =
            _mm256_add_ps(xi_46, _mm256_set_ps(xi_61, xi_61, xi_61, xi_61,
                                               xi_61, xi_61, xi_61, xi_61));
        const __m256d xi_71 = _mm256_add_ps(
            xi_67,
            _mm256_set_ps(
                _mm256_mul_ps(xi_65, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_65, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_65, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_65, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_65, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_65, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_65, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_65, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0))));
        const __m256d xi_72 = _mm256_add_ps(xi_64, xi_71);
        const __m256d xi_73 = _mm256_mul_ps(
            _mm256_set_ps(_mm256_add_ps(_mm256_mul_ps(u_2, xi_57),
                                        _mm256_mul_ps(xi_261, xi_56)),
                          _mm256_add_ps(_mm256_mul_ps(u_2, xi_57),
                                        _mm256_mul_ps(xi_261, xi_56)),
                          _mm256_add_ps(_mm256_mul_ps(u_2, xi_57),
                                        _mm256_mul_ps(xi_261, xi_56)),
                          _mm256_add_ps(_mm256_mul_ps(u_2, xi_57),
                                        _mm256_mul_ps(xi_261, xi_56)),
                          _mm256_add_ps(_mm256_mul_ps(u_2, xi_57),
                                        _mm256_mul_ps(xi_261, xi_56)),
                          _mm256_add_ps(_mm256_mul_ps(u_2, xi_57),
                                        _mm256_mul_ps(xi_261, xi_56)),
                          _mm256_add_ps(_mm256_mul_ps(u_2, xi_57),
                                        _mm256_mul_ps(xi_261, xi_56)),
                          _mm256_add_ps(_mm256_mul_ps(u_2, xi_57),
                                        _mm256_mul_ps(xi_261, xi_56))),
            _mm256_set_ps(xi_58, xi_58, xi_58, xi_58, xi_58, xi_58, xi_58,
                          xi_58));
        const __m256d xi_74 = _mm256_mul_ps(
            _mm256_mul_ps(
                _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0),
                _mm256_set_ps(xi_41, xi_41, xi_41, xi_41, xi_41, xi_41, xi_41,
                              xi_41)),
            _mm256_set_ps(xi_63, xi_63, xi_63, xi_63, xi_63, xi_63, xi_63,
                          xi_63));
        const __m256d xi_76 =
            _mm256_add_ps(xi_44, _mm256_set_ps(xi_75, xi_75, xi_75, xi_75,
                                               xi_75, xi_75, xi_75, xi_75));
        const __m256d xi_77 = _mm256_add_ps(xi_74, xi_76);
        const __m256d xi_78 =
            _mm256_mul_ps(xi_73, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                               -1.0, -1.0, -1.0));
        const __m256d xi_79 = _mm256_mul_ps(
            _mm256_add_ps(
                _mm256_mul_ps(
                    _mm256_mul_ps(
                        _mm256_set_ps(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                        _mm256_set_ps(u_0, u_0, u_0, u_0, u_0, u_0, u_0, u_0)),
                    _mm256_set_ps(xi_261, xi_261, xi_261, xi_261, xi_261,
                                  xi_261, xi_261, xi_261)),
                _mm256_mul_ps(
                    _mm256_mul_ps(
                        _mm256_set_ps(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                        _mm256_set_ps(u_2, u_2, u_2, u_2, u_2, u_2, u_2, u_2)),
                    _mm256_set_ps(xi_269, xi_269, xi_269, xi_269, xi_269,
                                  xi_269, xi_269, xi_269))),
            _mm256_set_ps(xi_58, xi_58, xi_58, xi_58, xi_58, xi_58, xi_58,
                          xi_58));
        const __m256d xi_80 =
            _mm256_mul_ps(xi_79, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                               -1.0, -1.0, -1.0));
        const __m256d xi_81 = _mm256_mul_ps(
            _mm256_mul_ps(
                _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0),
                _mm256_set_ps(xi_40, xi_40, xi_40, xi_40, xi_40, xi_40, xi_40,
                              xi_40)),
            _mm256_set_ps(xi_63, xi_63, xi_63, xi_63, xi_63, xi_63, xi_63,
                          xi_63));
        const __m256d xi_82 = _mm256_add_ps(_mm256_add_ps(xi_66, xi_76), xi_81);
        const __m256d xi_83 = _mm256_add_ps(
            xi_44,
            _mm256_set_ps(
                _mm256_mul_ps(xi_75, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_75, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_75, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_75, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_75, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_75, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_75, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_75, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0))));
        const __m256d xi_84 = _mm256_add_ps(xi_74, xi_83);
        const __m256d xi_85 = _mm256_add_ps(_mm256_add_ps(xi_66, xi_81), xi_83);
        const __m256 xi_100 = _mm256_mul_ps(rho, (_mm256_mul_ps(u_2, u_2)));
        const __m256d xi_108 = _mm256_mul_ps(
            _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk, omega_bulk,
                          omega_bulk, omega_bulk, omega_bulk, omega_bulk),
            _mm256_set_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_add_ps(_mm256_add_ps(xi_100, xi_102),
                                                  xi_105),
                                    xi_107),
                                xi_17),
                            xi_22),
                        xi_259),
                    xi_99),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_add_ps(_mm256_add_ps(xi_100, xi_102),
                                                  xi_105),
                                    xi_107),
                                xi_17),
                            xi_22),
                        xi_259),
                    xi_99),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_add_ps(_mm256_add_ps(xi_100, xi_102),
                                                  xi_105),
                                    xi_107),
                                xi_17),
                            xi_22),
                        xi_259),
                    xi_99),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_add_ps(_mm256_add_ps(xi_100, xi_102),
                                                  xi_105),
                                    xi_107),
                                xi_17),
                            xi_22),
                        xi_259),
                    xi_99),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_add_ps(_mm256_add_ps(xi_100, xi_102),
                                                  xi_105),
                                    xi_107),
                                xi_17),
                            xi_22),
                        xi_259),
                    xi_99),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_add_ps(_mm256_add_ps(xi_100, xi_102),
                                                  xi_105),
                                    xi_107),
                                xi_17),
                            xi_22),
                        xi_259),
                    xi_99),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_add_ps(_mm256_add_ps(xi_100, xi_102),
                                                  xi_105),
                                    xi_107),
                                xi_17),
                            xi_22),
                        xi_259),
                    xi_99),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_add_ps(_mm256_add_ps(xi_100, xi_102),
                                                  xi_105),
                                    xi_107),
                                xi_17),
                            xi_22),
                        xi_259),
                    xi_99)));
        const __m256 xi_143 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_mul_ps(xi_100, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                xi_253),
            xi_264);
        const __m256d xi_144 = _mm256_mul_ps(
            _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear,
                          omega_shear, omega_shear, omega_shear, omega_shear),
            _mm256_set_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_mul_ps(
                                        xi_262,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                                    xi_0),
                                xi_105),
                            xi_142),
                        xi_143),
                    xi_16),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_mul_ps(
                                        xi_262,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                                    xi_0),
                                xi_105),
                            xi_142),
                        xi_143),
                    xi_16),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_mul_ps(
                                        xi_262,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                                    xi_0),
                                xi_105),
                            xi_142),
                        xi_143),
                    xi_16),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_mul_ps(
                                        xi_262,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                                    xi_0),
                                xi_105),
                            xi_142),
                        xi_143),
                    xi_16),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_mul_ps(
                                        xi_262,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                                    xi_0),
                                xi_105),
                            xi_142),
                        xi_143),
                    xi_16),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_mul_ps(
                                        xi_262,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                                    xi_0),
                                xi_105),
                            xi_142),
                        xi_143),
                    xi_16),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_mul_ps(
                                        xi_262,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                                    xi_0),
                                xi_105),
                            xi_142),
                        xi_143),
                    xi_16),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_mul_ps(
                                        xi_262,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                                    xi_0),
                                xi_105),
                            xi_142),
                        xi_143),
                    xi_16)));
        const __m256d xi_145 =
            _mm256_mul_ps(xi_144, _mm256_set_ps(0.125, 0.125, 0.125, 0.125,
                                                0.125, 0.125, 0.125, 0.125));
        const __m256d xi_147 = _mm256_mul_ps(
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
                                                    _mm256_mul_ps(
                                                        _mm256_set_ps(
                                                            -2.0, -2.0, -2.0,
                                                            -2.0, -2.0, -2.0,
                                                            -2.0, -2.0),
                                                        _mm256_set_ps(
                                                            xi_257, xi_257,
                                                            xi_257, xi_257,
                                                            xi_257, xi_257,
                                                            xi_257, xi_257)),
                                                    _mm256_mul_ps(
                                                        _mm256_set_ps(
                                                            -2.0, -2.0, -2.0,
                                                            -2.0, -2.0, -2.0,
                                                            -2.0, -2.0),
                                                        _mm256_set_ps(
                                                            xi_263, xi_263,
                                                            xi_263, xi_263,
                                                            xi_263, xi_263,
                                                            xi_263, xi_263))),
                                                _mm256_mul_ps(
                                                    _mm256_set_ps(2.0, 2.0, 2.0,
                                                                  2.0, 2.0, 2.0,
                                                                  2.0, 2.0),
                                                    _mm256_set_ps(
                                                        xi_99, xi_99, xi_99,
                                                        xi_99, xi_99, xi_99,
                                                        xi_99, xi_99))),
                                            _mm256_set_ps(
                                                _mm256_mul_ps(
                                                    xi_104,
                                                    _mm256_set_ps(-1.0, -1.0,
                                                                  -1.0, -1.0,
                                                                  -1.0, -1.0,
                                                                  -1.0, -1.0)),
                                                _mm256_mul_ps(
                                                    xi_104,
                                                    _mm256_set_ps(-1.0, -1.0,
                                                                  -1.0, -1.0,
                                                                  -1.0, -1.0,
                                                                  -1.0, -1.0)),
                                                _mm256_mul_ps(
                                                    xi_104,
                                                    _mm256_set_ps(-1.0, -1.0,
                                                                  -1.0, -1.0,
                                                                  -1.0, -1.0,
                                                                  -1.0, -1.0)),
                                                _mm256_mul_ps(
                                                    xi_104,
                                                    _mm256_set_ps(-1.0, -1.0,
                                                                  -1.0, -1.0,
                                                                  -1.0, -1.0,
                                                                  -1.0, -1.0)),
                                                _mm256_mul_ps(
                                                    xi_104,
                                                    _mm256_set_ps(-1.0, -1.0,
                                                                  -1.0, -1.0,
                                                                  -1.0, -1.0,
                                                                  -1.0, -1.0)),
                                                _mm256_mul_ps(
                                                    xi_104,
                                                    _mm256_set_ps(-1.0, -1.0,
                                                                  -1.0, -1.0,
                                                                  -1.0, -1.0,
                                                                  -1.0, -1.0)),
                                                _mm256_mul_ps(
                                                    xi_104,
                                                    _mm256_set_ps(-1.0, -1.0,
                                                                  -1.0, -1.0,
                                                                  -1.0, -1.0,
                                                                  -1.0, -1.0)),
                                                _mm256_mul_ps(
                                                    xi_104,
                                                    _mm256_set_ps(
                                                        -1.0, -1.0, -1.0, -1.0,
                                                        -1.0, -1.0, -1.0,
                                                        -1.0)))),
                                        _mm256_set_ps(xi_103, xi_103, xi_103,
                                                      xi_103, xi_103, xi_103,
                                                      xi_103, xi_103)),
                                    _mm256_set_ps(xi_107, xi_107, xi_107,
                                                  xi_107, xi_107, xi_107,
                                                  xi_107, xi_107)),
                                _mm256_set_ps(xi_143, xi_143, xi_143, xi_143,
                                              xi_143, xi_143, xi_143, xi_143)),
                            _mm256_set_ps(xi_254, xi_254, xi_254, xi_254,
                                          xi_254, xi_254, xi_254, xi_254)),
                        _mm256_set_ps(xi_262, xi_262, xi_262, xi_262, xi_262,
                                      xi_262, xi_262, xi_262)),
                    _mm256_set_ps(xi_9, xi_9, xi_9, xi_9, xi_9, xi_9, xi_9,
                                  xi_9)),
                _mm256_set_ps(xi_96, xi_96, xi_96, xi_96, xi_96, xi_96, xi_96,
                              xi_96)),
            _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear,
                          omega_shear, omega_shear, omega_shear, omega_shear));
        const __m256d xi_149 = _mm256_add_ps(
            _mm256_mul_ps(
                xi_147,
                _mm256_set_ps(-0.0416666666666667, -0.0416666666666667,
                              -0.0416666666666667, -0.0416666666666667,
                              -0.0416666666666667, -0.0416666666666667,
                              -0.0416666666666667, -0.0416666666666667)),
            _mm256_mul_ps(
                xi_148, _mm256_set_ps(-0.166666666666667, -0.166666666666667,
                                      -0.166666666666667, -0.166666666666667,
                                      -0.166666666666667, -0.166666666666667,
                                      -0.166666666666667, -0.166666666666667)));
        const __m256d xi_150 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_mul_ps(xi_89, _mm256_set_ps(-0.1, -0.1, -0.1, -0.1, -0.1,
                                                   -0.1, -0.1, -0.1)),
                _mm256_mul_ps(xi_95,
                              _mm256_set_ps(-0.05, -0.05, -0.05, -0.05, -0.05,
                                            -0.05, -0.05, -0.05))),
            xi_149);
        const __m256d xi_151 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_mul_ps(
                                xi_92,
                                _mm256_set_ps(
                                    0.0285714285714286, 0.0285714285714286,
                                    0.0285714285714286, 0.0285714285714286,
                                    0.0285714285714286, 0.0285714285714286,
                                    0.0285714285714286, 0.0285714285714286)),
                            _mm256_mul_ps(
                                xi_98,
                                _mm256_set_ps(
                                    0.0142857142857143, 0.0142857142857143,
                                    0.0142857142857143, 0.0142857142857143,
                                    0.0142857142857143, 0.0142857142857143,
                                    0.0142857142857143, 0.0142857142857143))),
                        xi_141),
                    xi_145),
                xi_146),
            xi_150);
        const __m256d xi_166 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_mul_ps(
                            xi_147,
                            _mm256_set_ps(
                                0.0833333333333333, 0.0833333333333333,
                                0.0833333333333333, 0.0833333333333333,
                                0.0833333333333333, 0.0833333333333333,
                                0.0833333333333333, 0.0833333333333333)),
                        _mm256_mul_ps(
                            xi_148, _mm256_set_ps(
                                        0.333333333333333, 0.333333333333333,
                                        0.333333333333333, 0.333333333333333,
                                        0.333333333333333, 0.333333333333333,
                                        0.333333333333333, 0.333333333333333))),
                    _mm256_mul_ps(
                        xi_92, _mm256_set_ps(
                                   -0.0714285714285714, -0.0714285714285714,
                                   -0.0714285714285714, -0.0714285714285714,
                                   -0.0714285714285714, -0.0714285714285714,
                                   -0.0714285714285714, -0.0714285714285714))),
                _mm256_mul_ps(
                    xi_98,
                    _mm256_set_ps(-0.0357142857142857, -0.0357142857142857,
                                  -0.0357142857142857, -0.0357142857142857,
                                  -0.0357142857142857, -0.0357142857142857,
                                  -0.0357142857142857, -0.0357142857142857))),
            xi_146);
        const __m256 xi_171 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_mul_ps(rho, u_2),
                                _mm256_mul_ps(vel2Term,
                                              _mm256_set_ps(-1.0, -1.0, -1.0,
                                                            -1.0, -1.0, -1.0,
                                                            -1.0, -1.0))),
                            xi_101),
                        xi_106),
                    xi_168),
                xi_252),
            xi_6);
        const __m256 xi_172 = _mm256_mul_ps(
            xi_171, _mm256_set_ps(xi_120, xi_120, xi_120, xi_120, xi_120,
                                  xi_120, xi_120, xi_120));
        const __m256d xi_180 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_mul_ps(
                                    xi_112,
                                    _mm256_set_ps(
                                        0.0158730158730159, 0.0158730158730159,
                                        0.0158730158730159, 0.0158730158730159,
                                        0.0158730158730159, 0.0158730158730159,
                                        0.0158730158730159,
                                        0.0158730158730159)),
                                _mm256_mul_ps(xi_141,
                                              _mm256_set_ps(-1.0, -1.0, -1.0,
                                                            -1.0, -1.0, -1.0,
                                                            -1.0, -1.0))),
                            _mm256_mul_ps(
                                xi_145, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0))),
                        _mm256_mul_ps(
                            xi_88,
                            _mm256_set_ps(
                                0.0952380952380952, 0.0952380952380952,
                                0.0952380952380952, 0.0952380952380952,
                                0.0952380952380952, 0.0952380952380952,
                                0.0952380952380952, 0.0952380952380952))),
                    _mm256_mul_ps(
                        xi_92, _mm256_set_ps(
                                   -0.0428571428571429, -0.0428571428571429,
                                   -0.0428571428571429, -0.0428571428571429,
                                   -0.0428571428571429, -0.0428571428571429,
                                   -0.0428571428571429, -0.0428571428571429))),
                _mm256_mul_ps(
                    xi_98,
                    _mm256_set_ps(-0.0214285714285714, -0.0214285714285714,
                                  -0.0214285714285714, -0.0214285714285714,
                                  -0.0214285714285714, -0.0214285714285714,
                                  -0.0214285714285714, -0.0214285714285714))),
            xi_150);
        const __m256d xi_183 = _mm256_mul_ps(
            xi_144, _mm256_set_ps(0.0625, 0.0625, 0.0625, 0.0625, 0.0625,
                                  0.0625, 0.0625, 0.0625));
        const __m256d xi_188 = _mm256_add_ps(
            _mm256_mul_ps(
                xi_108, _mm256_set_ps(0.0416666666666667, 0.0416666666666667,
                                      0.0416666666666667, 0.0416666666666667,
                                      0.0416666666666667, 0.0416666666666667,
                                      0.0416666666666667, 0.0416666666666667)),
            _mm256_mul_ps(
                xi_91, _mm256_set_ps(0.0833333333333333, 0.0833333333333333,
                                     0.0833333333333333, 0.0833333333333333,
                                     0.0833333333333333, 0.0833333333333333,
                                     0.0833333333333333, 0.0833333333333333)));
        const __m256d xi_189 = _mm256_add_ps(
            xi_188, _mm256_set_ps(xi_187, xi_187, xi_187, xi_187, xi_187,
                                  xi_187, xi_187, xi_187));
        const __m256d xi_190 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(xi_152, xi_182), xi_183),
                    xi_184),
                xi_185),
            xi_189);
        const __m256d xi_192 = _mm256_add_ps(
            _mm256_mul_ps(
                xi_147, _mm256_set_ps(0.0208333333333333, 0.0208333333333333,
                                      0.0208333333333333, 0.0208333333333333,
                                      0.0208333333333333, 0.0208333333333333,
                                      0.0208333333333333, 0.0208333333333333)),
            _mm256_mul_ps(
                xi_148, _mm256_set_ps(0.0833333333333333, 0.0833333333333333,
                                      0.0833333333333333, 0.0833333333333333,
                                      0.0833333333333333, 0.0833333333333333,
                                      0.0833333333333333, 0.0833333333333333)));
        const __m256d xi_193 = _mm256_add_ps(
            xi_192,
            _mm256_set_ps(
                _mm256_mul_ps(xi_191, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_191, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_191, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_191, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_191, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_191, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_191, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_191, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0))));
        const __m256d xi_194 = _mm256_add_ps(xi_167, xi_193);
        const __m256d xi_200 = _mm256_add_ps(
            xi_192, _mm256_set_ps(xi_191, xi_191, xi_191, xi_191, xi_191,
                                  xi_191, xi_191, xi_191));
        const __m256d xi_201 = _mm256_add_ps(xi_165, xi_200);
        const __m256d xi_202 = _mm256_add_ps(
            xi_188,
            _mm256_set_ps(
                _mm256_mul_ps(xi_187, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_187, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_187, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_187, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_187, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_187, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_187, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_187, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0))));
        const __m256d xi_203 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(xi_137, xi_182), xi_183),
                    xi_184),
                xi_185),
            xi_202);
        const __m256d xi_206 = _mm256_mul_ps(
            _mm256_set_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_mul_ps(u_2, xi_117), xi_113),
                        xi_17),
                    xi_251),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_mul_ps(u_2, xi_117), xi_113),
                        xi_17),
                    xi_251),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_mul_ps(u_2, xi_117), xi_113),
                        xi_17),
                    xi_251),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_mul_ps(u_2, xi_117), xi_113),
                        xi_17),
                    xi_251),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_mul_ps(u_2, xi_117), xi_113),
                        xi_17),
                    xi_251),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_mul_ps(u_2, xi_117), xi_113),
                        xi_17),
                    xi_251),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_mul_ps(u_2, xi_117), xi_113),
                        xi_17),
                    xi_251),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_mul_ps(u_2, xi_117), xi_113),
                        xi_17),
                    xi_251)),
            _mm256_set_ps(xi_196, xi_196, xi_196, xi_196, xi_196, xi_196,
                          xi_196, xi_196));
        const __m256d xi_210 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_add_ps(xi_149, xi_204), xi_205),
                        xi_206),
                    xi_207),
                xi_208),
            xi_209);
        const __m256 xi_219 = _mm256_mul_ps(
            xi_171, _mm256_set_ps(xi_186, xi_186, xi_186, xi_186, xi_186,
                                  xi_186, xi_186, xi_186));
        const __m256d xi_221 = _mm256_add_ps(
            xi_220, _mm256_set_ps(xi_219, xi_219, xi_219, xi_219, xi_219,
                                  xi_219, xi_219, xi_219));
        const __m256d xi_222 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_mul_ps(xi_212,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_216,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0))),
                    xi_214),
                xi_218),
            xi_221);
        const __m256d xi_227 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_mul_ps(xi_223,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_225,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0))),
                    xi_189),
                xi_224),
            xi_226);
        const __m256d xi_228 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_mul_ps(xi_224,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_226,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0))),
                    xi_202),
                xi_223),
            xi_225);
        const __m256d xi_229 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_mul_ps(xi_204,
                                              _mm256_set_ps(-1.0, -1.0, -1.0,
                                                            -1.0, -1.0, -1.0,
                                                            -1.0, -1.0)),
                                _mm256_mul_ps(xi_206,
                                              _mm256_set_ps(-1.0, -1.0, -1.0,
                                                            -1.0, -1.0, -1.0,
                                                            -1.0, -1.0))),
                            xi_149),
                        xi_205),
                    xi_207),
                xi_208),
            xi_209);
        const __m256d xi_231 =
            _mm256_mul_ps(xi_183, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256d xi_234 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_add_ps(xi_181, xi_188), xi_221),
                        xi_230),
                    xi_231),
                xi_232),
            xi_233);
        const __m256d xi_236 = _mm256_mul_ps(
            _mm256_set_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_mul_ps(u_2, xi_153), xi_10),
                        xi_156),
                    xi_252),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_mul_ps(u_2, xi_153), xi_10),
                        xi_156),
                    xi_252),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_mul_ps(u_2, xi_153), xi_10),
                        xi_156),
                    xi_252),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_mul_ps(u_2, xi_153), xi_10),
                        xi_156),
                    xi_252),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_mul_ps(u_2, xi_153), xi_10),
                        xi_156),
                    xi_252),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_mul_ps(u_2, xi_153), xi_10),
                        xi_156),
                    xi_252),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_mul_ps(u_2, xi_153), xi_10),
                        xi_156),
                    xi_252),
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_mul_ps(u_2, xi_153), xi_10),
                        xi_156),
                    xi_252)),
            _mm256_set_ps(xi_196, xi_196, xi_196, xi_196, xi_196, xi_196,
                          xi_196, xi_196));
        const __m256d xi_237 = _mm256_add_ps(
            _mm256_mul_ps(xi_235, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)),
            _mm256_mul_ps(xi_236, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)));
        const __m256d xi_242 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_mul_ps(xi_238,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_240,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0))),
                    xi_193),
                xi_239),
            xi_241);
        const __m256d xi_243 = _mm256_add_ps(xi_235, xi_236);
        const __m256d xi_244 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_mul_ps(xi_239,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_241,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0))),
                    xi_200),
                xi_238),
            xi_240);
        const __m256d xi_245 = _mm256_add_ps(
            xi_220,
            _mm256_set_ps(
                _mm256_mul_ps(xi_219, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_219, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_219, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_219, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_219, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_219, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_219, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_219, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0))));
        const __m256d xi_246 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_mul_ps(xi_214,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_218,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0))),
                    xi_212),
                xi_216),
            xi_245);
        const __m256d xi_247 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_add_ps(xi_179, xi_188), xi_230),
                        xi_231),
                    xi_232),
                xi_233),
            xi_245);
        const __m256d forceTerm_0 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_mul_ps(xi_32, _mm256_set_ps(-1.5, -1.5, -1.5, -1.5,
                                                       -1.5, -1.5, -1.5, -1.5)),
                    _mm256_mul_ps(
                        _mm256_mul_ps(_mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0),
                                      _mm256_set_ps(xi_35, xi_35, xi_35, xi_35,
                                                    xi_35, xi_35, xi_35,
                                                    xi_35)),
                        _mm256_set_ps(xi_38, xi_38, xi_38, xi_38, xi_38, xi_38,
                                      xi_38, xi_38))),
                _mm256_mul_ps(
                    _mm256_mul_ps(_mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0),
                                  _mm256_set_ps(xi_38, xi_38, xi_38, xi_38,
                                                xi_38, xi_38, xi_38, xi_38)),
                    _mm256_set_ps(xi_40, xi_40, xi_40, xi_40, xi_40, xi_40,
                                  xi_40, xi_40))),
            _mm256_mul_ps(
                _mm256_mul_ps(_mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
                                            -1.0, -1.0),
                              _mm256_set_ps(xi_38, xi_38, xi_38, xi_38, xi_38,
                                            xi_38, xi_38, xi_38)),
                _mm256_set_ps(xi_41, xi_41, xi_41, xi_41, xi_41, xi_41, xi_41,
                              xi_41)));
        const __m256d forceTerm_1 =
            _mm256_add_ps(xi_49, _mm256_set_ps(xi_42, xi_42, xi_42, xi_42,
                                               xi_42, xi_42, xi_42, xi_42));
        const __m256d forceTerm_2 = _mm256_add_ps(
            xi_49,
            _mm256_set_ps(
                _mm256_mul_ps(xi_42, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_42, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_42, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_42, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_42, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_42, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_42, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_42, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0))));
        const __m256d forceTerm_3 = _mm256_add_ps(
            xi_53,
            _mm256_set_ps(
                _mm256_mul_ps(xi_50, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_50, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_50, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_50, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_50, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_50, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_50, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_50, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0))));
        const __m256d forceTerm_4 =
            _mm256_add_ps(xi_53, _mm256_set_ps(xi_50, xi_50, xi_50, xi_50,
                                               xi_50, xi_50, xi_50, xi_50));
        const __m256d forceTerm_5 =
            _mm256_add_ps(xi_55, _mm256_set_ps(xi_54, xi_54, xi_54, xi_54,
                                               xi_54, xi_54, xi_54, xi_54));
        const __m256d forceTerm_6 = _mm256_add_ps(
            xi_55,
            _mm256_set_ps(
                _mm256_mul_ps(xi_54, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_54, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_54, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_54, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_54, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_54, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_54, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_mul_ps(xi_54, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0))));
        const __m256d forceTerm_7 =
            _mm256_add_ps(_mm256_add_ps(xi_60, xi_62), xi_69);
        const __m256d forceTerm_8 =
            _mm256_add_ps(_mm256_add_ps(xi_59, xi_69), xi_70);
        const __m256d forceTerm_9 =
            _mm256_add_ps(_mm256_add_ps(xi_59, xi_62), xi_72);
        const __m256d forceTerm_10 =
            _mm256_add_ps(_mm256_add_ps(xi_60, xi_70), xi_72);
        const __m256d forceTerm_11 =
            _mm256_add_ps(_mm256_add_ps(xi_68, xi_73), xi_77);
        const __m256d forceTerm_12 =
            _mm256_add_ps(_mm256_add_ps(xi_71, xi_77), xi_78);
        const __m256d forceTerm_13 =
            _mm256_add_ps(_mm256_add_ps(xi_62, xi_80), xi_82);
        const __m256d forceTerm_14 =
            _mm256_add_ps(_mm256_add_ps(xi_70, xi_79), xi_82);
        const __m256d forceTerm_15 =
            _mm256_add_ps(_mm256_add_ps(xi_68, xi_78), xi_84);
        const __m256d forceTerm_16 =
            _mm256_add_ps(_mm256_add_ps(xi_71, xi_73), xi_84);
        const __m256d forceTerm_17 =
            _mm256_add_ps(_mm256_add_ps(xi_62, xi_79), xi_85);
        const __m256d forceTerm_18 =
            _mm256_add_ps(_mm256_add_ps(xi_70, xi_80), xi_85);
        _mm256_store_ps(
            &_data_pdfs_20_30_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_add_ps(
                                        _mm256_add_ps(
                                            _mm256_add_ps(
                                                _mm256_mul_ps(
                                                    xi_108,
                                                    _mm256_set_ps(-0.5, -0.5,
                                                                  -0.5, -0.5,
                                                                  -0.5, -0.5,
                                                                  -0.5, -0.5)),
                                                _mm256_mul_ps(
                                                    xi_112,
                                                    _mm256_set_ps(
                                                        0.0238095238095238,
                                                        0.0238095238095238,
                                                        0.0238095238095238,
                                                        0.0238095238095238,
                                                        0.0238095238095238,
                                                        0.0238095238095238,
                                                        0.0238095238095238,
                                                        0.0238095238095238))),
                                            _mm256_mul_ps(
                                                xi_88, _mm256_set_ps(
                                                           0.142857142857143,
                                                           0.142857142857143,
                                                           0.142857142857143,
                                                           0.142857142857143,
                                                           0.142857142857143,
                                                           0.142857142857143,
                                                           0.142857142857143,
                                                           0.142857142857143))),
                                        _mm256_mul_ps(
                                            xi_89,
                                            _mm256_set_ps(0.2, 0.2, 0.2, 0.2,
                                                          0.2, 0.2, 0.2, 0.2))),
                                    _mm256_mul_ps(
                                        xi_91,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0))),
                                _mm256_mul_ps(
                                    xi_92,
                                    _mm256_set_ps(
                                        0.0857142857142857, 0.0857142857142857,
                                        0.0857142857142857, 0.0857142857142857,
                                        0.0857142857142857, 0.0857142857142857,
                                        0.0857142857142857,
                                        0.0857142857142857))),
                            _mm256_mul_ps(xi_95,
                                          _mm256_set_ps(0.1, 0.1, 0.1, 0.1, 0.1,
                                                        0.1, 0.1, 0.1))),
                        _mm256_mul_ps(
                            xi_98,
                            _mm256_set_ps(
                                0.0428571428571429, 0.0428571428571429,
                                0.0428571428571429, 0.0428571428571429,
                                0.0428571428571429, 0.0428571428571429,
                                0.0428571428571429, 0.0428571428571429))),
                    forceTerm_0),
                _mm256_set_ps(xi_259, xi_259, xi_259, xi_259, xi_259, xi_259,
                              xi_259, xi_259)));
        _mm256_store_ps(
            &_data_pdfs_20_31_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_mul_ps(
                                        xi_116,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                                    _mm256_mul_ps(
                                        xi_126,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0))),
                                forceTerm_1),
                            xi_137),
                        xi_151),
                    _mm256_set_ps(xi_121, xi_121, xi_121, xi_121, xi_121,
                                  xi_121, xi_121, xi_121)),
                _mm256_set_ps(xi_262, xi_262, xi_262, xi_262, xi_262, xi_262,
                              xi_262, xi_262)));
        _mm256_store_ps(
            &_data_pdfs_20_32_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(_mm256_add_ps(forceTerm_2, xi_116),
                                          xi_126),
                            xi_151),
                        xi_152),
                    _mm256_set_ps(
                        _mm256_mul_ps(xi_121,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_121,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_121,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_121,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_121,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_121,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_121,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_121,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)))),
                _mm256_set_ps(xi_254, xi_254, xi_254, xi_254, xi_254, xi_254,
                              xi_254, xi_254)));
        _mm256_store_ps(
            &_data_pdfs_20_33_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(_mm256_add_ps(forceTerm_3, xi_158),
                                          xi_160),
                            xi_165),
                        xi_166),
                    _mm256_set_ps(
                        _mm256_mul_ps(xi_155,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_155,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_155,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_155,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_155,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_155,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_155,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_155,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)))),
                _mm256_set_ps(xi_263, xi_263, xi_263, xi_263, xi_263, xi_263,
                              xi_263, xi_263)));
        _mm256_store_ps(
            &_data_pdfs_20_34_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_mul_ps(
                                        xi_158,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                                    _mm256_mul_ps(
                                        xi_160,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0))),
                                forceTerm_4),
                            xi_166),
                        xi_167),
                    _mm256_set_ps(xi_155, xi_155, xi_155, xi_155, xi_155,
                                  xi_155, xi_155, xi_155)),
                _mm256_set_ps(xi_257, xi_257, xi_257, xi_257, xi_257, xi_257,
                              xi_257, xi_257)));
        _mm256_store_ps(
            &_data_pdfs_20_35_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_mul_ps(
                                        xi_170,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                                    _mm256_mul_ps(
                                        xi_174,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0))),
                                forceTerm_5),
                            xi_179),
                        xi_180),
                    _mm256_set_ps(xi_172, xi_172, xi_172, xi_172, xi_172,
                                  xi_172, xi_172, xi_172)),
                _mm256_set_ps(xi_253, xi_253, xi_253, xi_253, xi_253, xi_253,
                              xi_253, xi_253)));
        _mm256_store_ps(
            &_data_pdfs_20_36_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(_mm256_add_ps(forceTerm_6, xi_170),
                                          xi_174),
                            xi_180),
                        xi_181),
                    _mm256_set_ps(
                        _mm256_mul_ps(xi_172,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_172,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_172,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_172,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_172,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_172,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_172,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_172,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)))),
                _mm256_set_ps(xi_264, xi_264, xi_264, xi_264, xi_264, xi_264,
                              xi_264, xi_264)));
        _mm256_store_ps(
            &_data_pdfs_20_37_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(forceTerm_7, xi_190), xi_194),
                    xi_198),
                _mm256_set_ps(xi_249, xi_249, xi_249, xi_249, xi_249, xi_249,
                              xi_249, xi_249)));
        _mm256_store_ps(
            &_data_pdfs_20_38_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(forceTerm_8, xi_190), xi_199),
                    xi_201),
                _mm256_set_ps(xi_258, xi_258, xi_258, xi_258, xi_258, xi_258,
                              xi_258, xi_258)));
        _mm256_store_ps(
            &_data_pdfs_20_39_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(forceTerm_9, xi_194), xi_199),
                    xi_203),
                _mm256_set_ps(xi_268, xi_268, xi_268, xi_268, xi_268, xi_268,
                              xi_268, xi_268)));
        _mm256_store_ps(
            &_data_pdfs_20_310_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(forceTerm_10, xi_198), xi_201),
                    xi_203),
                _mm256_set_ps(xi_256, xi_256, xi_256, xi_256, xi_256, xi_256,
                              xi_256, xi_256)));
        _mm256_store_ps(
            &_data_pdfs_20_311_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(forceTerm_11, xi_210), xi_222),
                    xi_227),
                _mm256_set_ps(xi_265, xi_265, xi_265, xi_265, xi_265, xi_265,
                              xi_265, xi_265)));
        _mm256_store_ps(
            &_data_pdfs_20_312_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(forceTerm_12, xi_222), xi_228),
                    xi_229),
                _mm256_set_ps(xi_250, xi_250, xi_250, xi_250, xi_250, xi_250,
                              xi_250, xi_250)));
        _mm256_store_ps(
            &_data_pdfs_20_313_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(forceTerm_13, xi_234), xi_237),
                    xi_242),
                _mm256_set_ps(xi_255, xi_255, xi_255, xi_255, xi_255, xi_255,
                              xi_255, xi_255)));
        _mm256_store_ps(
            &_data_pdfs_20_314_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(forceTerm_14, xi_234), xi_243),
                    xi_244),
                _mm256_set_ps(xi_266, xi_266, xi_266, xi_266, xi_266, xi_266,
                              xi_266, xi_266)));
        _mm256_store_ps(
            &_data_pdfs_20_315_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(forceTerm_15, xi_227), xi_229),
                    xi_246),
                _mm256_set_ps(xi_251, xi_251, xi_251, xi_251, xi_251, xi_251,
                              xi_251, xi_251)));
        _mm256_store_ps(
            &_data_pdfs_20_316_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(forceTerm_16, xi_210), xi_228),
                    xi_246),
                _mm256_set_ps(xi_267, xi_267, xi_267, xi_267, xi_267, xi_267,
                              xi_267, xi_267)));
        _mm256_store_ps(
            &_data_pdfs_20_317_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(forceTerm_17, xi_242), xi_243),
                    xi_247),
                _mm256_set_ps(xi_260, xi_260, xi_260, xi_260, xi_260, xi_260,
                              xi_260, xi_260)));
        _mm256_store_ps(
            &_data_pdfs_20_318_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(forceTerm_18, xi_237), xi_244),
                    xi_247),
                _mm256_set_ps(xi_252, xi_252, xi_252, xi_252, xi_252, xi_252,
                              xi_252, xi_252)));
      }
    }
  }
}
} // namespace internal_collidesweepsingleprecisionthermalizedavx

void CollideSweepSinglePrecisionThermalizedAVX::operator()(IBlock *block) {
  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);
  auto force = block->getData<field::GhostLayerField<float, 3>>(forceID);

  auto &time_step = this->time_step_;
  auto &omega_odd = this->omega_odd_;
  auto block_offset_2 = this->block_offset_2_;
  auto &omega_shear = this->omega_shear_;
  auto block_offset_1 = this->block_offset_1_;
  auto block_offset_0 = this->block_offset_0_;
  auto &seed = this->seed_;
  auto &kT = this->kT_;
  auto &omega_even = this->omega_even_;
  auto &omega_bulk = this->omega_bulk_;
  block_offset_generator(block, block_offset_0, block_offset_1, block_offset_2);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force->nrOfGhostLayers()));
  float *RESTRICT const _data_force = force->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()));
  float *RESTRICT _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
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
  internal_collidesweepsingleprecisionthermalizedavx::
      collidesweepsingleprecisionthermalizedavx(
          _data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2,
          _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_1,
          _stride_pdfs_2, _stride_pdfs_3, block_offset_0, block_offset_1,
          block_offset_2, kT, omega_bulk, omega_even, omega_odd, omega_shear,
          seed, time_step);
}

void CollideSweepSinglePrecisionThermalizedAVX::runOnCellInterval(
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

  auto &time_step = this->time_step_;
  auto &omega_odd = this->omega_odd_;
  auto block_offset_2 = this->block_offset_2_;
  auto &omega_shear = this->omega_shear_;
  auto block_offset_1 = this->block_offset_1_;
  auto block_offset_0 = this->block_offset_0_;
  auto &seed = this->seed_;
  auto &kT = this->kT_;
  auto &omega_even = this->omega_even_;
  auto &omega_bulk = this->omega_bulk_;
  block_offset_generator(block, block_offset_0, block_offset_1, block_offset_2);
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
  internal_collidesweepsingleprecisionthermalizedavx::
      collidesweepsingleprecisionthermalizedavx(
          _data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2,
          _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_1,
          _stride_pdfs_2, _stride_pdfs_3, block_offset_0, block_offset_1,
          block_offset_2, kT, omega_bulk, omega_even, omega_odd, omega_shear,
          seed, time_step);
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