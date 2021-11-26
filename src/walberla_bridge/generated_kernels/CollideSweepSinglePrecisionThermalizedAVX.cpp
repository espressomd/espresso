// kernel generated with pystencils v0.4.3, lbmpy v0.4.3,
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

namespace internal_48c9ee502281a70505dce0378c55abd5 {
static FUNC_PREFIX void
collidesweepsingleprecisionthermalizedavx_collidesweepsingleprecisionthermalizedavx(
    float *RESTRICT const _data_force, float *RESTRICT _data_pdfs,
    int64_t const _size_force_0, int64_t const _size_force_1,
    int64_t const _size_force_2, int64_t const _stride_force_1,
    int64_t const _stride_force_2, int64_t const _stride_force_3,
    int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
    int64_t const _stride_pdfs_3, uint32_t block_offset_0,
    uint32_t block_offset_1, uint32_t block_offset_2, float kT,
    float omega_bulk, float omega_even, float omega_odd, float omega_shear,
    uint32_t seed, uint32_t time_step) {
  const float xi_29 = omega_bulk * 0.5f;
  const float xi_56 = omega_shear * 0.0416666666666667f;
  const float xi_61 = omega_bulk * 0.0416666666666667f;
  const float xi_72 = omega_shear * 0.125f;
  const float xi_110 = 2.4494897427831779f;
  const float xi_135 = omega_odd * 0.25f;
  const float xi_151 = omega_odd * 0.0833333333333333f;
  const float xi_216 = omega_shear * 0.25f;
  const float xi_231 = omega_odd * 0.0416666666666667f;
  const float xi_233 = omega_odd * 0.125f;
  const float rr_0 = 0.0f;
  const float xi_54 = rr_0 * 0.0416666666666667f;
  const float xi_140 = rr_0 * 0.166666666666667f;
  const float xi_206 = rr_0 * 0.0833333333333333f;
  for (int64_t ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1) {
    float *RESTRICT _data_pdfs_20_315 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_32 =
        _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    float *RESTRICT _data_pdfs_20_36 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_31 =
        _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    float *RESTRICT _data_pdfs_20_35 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_33 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    float *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    float *RESTRICT _data_pdfs_20_32 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_311 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_318 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_34 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_314 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_312 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_316 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_317 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_31 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_313 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1) {
      float *RESTRICT _data_pdfs_20_315_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_315;
      float *RESTRICT _data_force_20_32_10 =
          _stride_force_1 * ctr_1 + _data_force_20_32;
      float *RESTRICT _data_pdfs_20_36_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_36;
      float *RESTRICT _data_force_20_31_10 =
          _stride_force_1 * ctr_1 + _data_force_20_31;
      float *RESTRICT _data_pdfs_20_35_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_35;
      float *RESTRICT _data_pdfs_20_33_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      float *RESTRICT _data_pdfs_20_30_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      float *RESTRICT _data_pdfs_20_38_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_38;
      float *RESTRICT _data_force_20_30_10 =
          _stride_force_1 * ctr_1 + _data_force_20_30;
      float *RESTRICT _data_pdfs_20_32_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_32;
      float *RESTRICT _data_pdfs_20_37_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_37;
      float *RESTRICT _data_pdfs_20_39_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_39;
      float *RESTRICT _data_pdfs_20_311_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_311;
      float *RESTRICT _data_pdfs_20_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_318;
      float *RESTRICT _data_pdfs_20_34_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      float *RESTRICT _data_pdfs_20_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_314;
      float *RESTRICT _data_pdfs_20_312_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_312;
      float *RESTRICT _data_pdfs_20_310_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_310;
      float *RESTRICT _data_pdfs_20_316_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_316;
      float *RESTRICT _data_pdfs_20_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_317;
      float *RESTRICT _data_pdfs_20_31_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_31;
      float *RESTRICT _data_pdfs_20_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_313;
      for (int64_t ctr_0 = 0;
           ctr_0 < ((_size_force_0) % (8) == 0
                        ? _size_force_0
                        : ((int64_t)((_size_force_0) / (8)) + 1) * (8));
           ctr_0 += 8) {
        const __m256 xi_268 = _mm256_load_ps(&_data_pdfs_20_315_10[ctr_0]);
        const __m256 xi_269 = _mm256_load_ps(&_data_force_20_32_10[ctr_0]);
        const __m256 xi_270 = _mm256_load_ps(&_data_pdfs_20_36_10[ctr_0]);
        const __m256 xi_271 = _mm256_load_ps(&_data_force_20_31_10[ctr_0]);
        const __m256 xi_272 = _mm256_load_ps(&_data_pdfs_20_35_10[ctr_0]);
        const __m256 xi_273 = _mm256_load_ps(&_data_pdfs_20_33_10[ctr_0]);
        const __m256 xi_274 = _mm256_load_ps(&_data_pdfs_20_30_10[ctr_0]);
        const __m256 xi_275 = _mm256_load_ps(&_data_pdfs_20_38_10[ctr_0]);
        const __m256 xi_276 = _mm256_load_ps(&_data_force_20_30_10[ctr_0]);
        const __m256 xi_277 = _mm256_load_ps(&_data_pdfs_20_32_10[ctr_0]);
        const __m256 xi_278 = _mm256_load_ps(&_data_pdfs_20_37_10[ctr_0]);
        const __m256 xi_279 = _mm256_load_ps(&_data_pdfs_20_39_10[ctr_0]);
        const __m256 xi_280 = _mm256_load_ps(&_data_pdfs_20_311_10[ctr_0]);
        const __m256 xi_281 = _mm256_load_ps(&_data_pdfs_20_318_10[ctr_0]);
        const __m256 xi_282 = _mm256_load_ps(&_data_pdfs_20_34_10[ctr_0]);
        const __m256 xi_283 = _mm256_load_ps(&_data_pdfs_20_314_10[ctr_0]);
        const __m256 xi_284 = _mm256_load_ps(&_data_pdfs_20_312_10[ctr_0]);
        const __m256 xi_285 = _mm256_load_ps(&_data_pdfs_20_310_10[ctr_0]);
        const __m256 xi_286 = _mm256_load_ps(&_data_pdfs_20_316_10[ctr_0]);
        const __m256 xi_287 = _mm256_load_ps(&_data_pdfs_20_317_10[ctr_0]);
        const __m256 xi_288 = _mm256_load_ps(&_data_pdfs_20_31_10[ctr_0]);
        const __m256 xi_289 = _mm256_load_ps(&_data_pdfs_20_313_10[ctr_0]);

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

        const __m256 xi_0 = _mm256_add_ps(xi_281, xi_283);
        const __m256 xi_1 = _mm256_add_ps(xi_0, xi_282);
        const __m256 xi_2 =
            _mm256_add_ps(_mm256_add_ps(xi_268, xi_280), xi_288);
        const __m256 xi_3 = _mm256_add_ps(xi_272, xi_284);
        const __m256 xi_4 = _mm256_add_ps(xi_273, xi_279);
        const __m256 xi_5 = _mm256_add_ps(xi_277, xi_286);
        const __m256 xi_6 = _mm256_add_ps(xi_270, xi_287);
        const __m256 xi_9 =
            _mm256_mul_ps(xi_279, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_10 = _mm256_add_ps(
            _mm256_mul_ps(xi_278, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)),
            xi_9);
        const __m256 xi_11 =
            _mm256_mul_ps(xi_287, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_12 =
            _mm256_mul_ps(xi_289, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_13 =
            _mm256_mul_ps(xi_273, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_14 = _mm256_add_ps(_mm256_add_ps(xi_11, xi_12), xi_13);
        const __m256 xi_15 =
            _mm256_mul_ps(xi_277, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_16 =
            _mm256_mul_ps(xi_285, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_17 = _mm256_add_ps(xi_15, xi_16);
        const __m256 xi_18 =
            _mm256_mul_ps(xi_286, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_19 =
            _mm256_mul_ps(xi_284, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_20 = _mm256_add_ps(xi_18, xi_19);
        const __m256 xi_21 =
            _mm256_mul_ps(xi_281, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_22 = _mm256_add_ps(xi_11, xi_21);
        const __m256 xi_23 =
            _mm256_mul_ps(xi_268, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_24 =
            _mm256_mul_ps(xi_270, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_25 = _mm256_add_ps(
            _mm256_add_ps(_mm256_add_ps(xi_18, xi_23), xi_24), xi_280);
        const __m256 xi_30 = _mm256_mul_ps(
            xi_271, _mm256_set_ps(0.166666666666667f, 0.166666666666667f,
                                  0.166666666666667f, 0.166666666666667f,
                                  0.166666666666667f, 0.166666666666667f,
                                  0.166666666666667f, 0.166666666666667f));
        const __m256 xi_31 = _mm256_mul_ps(
            xi_271, _mm256_set_ps(0.0833333333333333f, 0.0833333333333333f,
                                  0.0833333333333333f, 0.0833333333333333f,
                                  0.0833333333333333f, 0.0833333333333333f,
                                  0.0833333333333333f, 0.0833333333333333f));
        const __m256 xi_43 = _mm256_mul_ps(
            xi_276, _mm256_set_ps(0.166666666666667f, 0.166666666666667f,
                                  0.166666666666667f, 0.166666666666667f,
                                  0.166666666666667f, 0.166666666666667f,
                                  0.166666666666667f, 0.166666666666667f));
        const __m256 xi_44 = _mm256_mul_ps(
            xi_276, _mm256_set_ps(0.0833333333333333f, 0.0833333333333333f,
                                  0.0833333333333333f, 0.0833333333333333f,
                                  0.0833333333333333f, 0.0833333333333333f,
                                  0.0833333333333333f, 0.0833333333333333f));
        const __m256 xi_50 = _mm256_mul_ps(
            xi_269, _mm256_set_ps(0.166666666666667f, 0.166666666666667f,
                                  0.166666666666667f, 0.166666666666667f,
                                  0.166666666666667f, 0.166666666666667f,
                                  0.166666666666667f, 0.166666666666667f));
        const __m256 xi_51 = _mm256_mul_ps(
            xi_269, _mm256_set_ps(0.0833333333333333f, 0.0833333333333333f,
                                  0.0833333333333333f, 0.0833333333333333f,
                                  0.0833333333333333f, 0.0833333333333333f,
                                  0.0833333333333333f, 0.0833333333333333f));
        const __m256 xi_68 =
            _mm256_mul_ps(xi_271, _mm256_set_ps(0.25f, 0.25f, 0.25f, 0.25f,
                                                0.25f, 0.25f, 0.25f, 0.25f));
        const __m256 xi_73 =
            _mm256_mul_ps(xi_271, _mm256_set_ps(xi_72, xi_72, xi_72, xi_72,
                                                xi_72, xi_72, xi_72, xi_72));
        const __m256 xi_113 =
            _mm256_mul_ps(xi_274, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_114 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_mul_ps(xi_270, _mm256_set_ps(3.0f, 3.0f, 3.0f, 3.0f,
                                                    3.0f, 3.0f, 3.0f, 3.0f)),
                _mm256_mul_ps(xi_272, _mm256_set_ps(3.0f, 3.0f, 3.0f, 3.0f,
                                                    3.0f, 3.0f, 3.0f, 3.0f))),
            xi_113);
        const __m256 xi_115 = _mm256_mul_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_mul_ps(xi_268,
                                                  _mm256_set_ps(-3.0f, -3.0f,
                                                                -3.0f, -3.0f,
                                                                -3.0f, -3.0f,
                                                                -3.0f, -3.0f)),
                                    _mm256_mul_ps(
                                        xi_277,
                                        _mm256_set_ps(3.0f, 3.0f, 3.0f, 3.0f,
                                                      3.0f, 3.0f, 3.0f, 3.0f))),
                                _mm256_mul_ps(xi_280,
                                              _mm256_set_ps(-3.0f, -3.0f, -3.0f,
                                                            -3.0f, -3.0f, -3.0f,
                                                            -3.0f, -3.0f))),
                            _mm256_mul_ps(xi_284,
                                          _mm256_set_ps(-3.0f, -3.0f, -3.0f,
                                                        -3.0f, -3.0f, -3.0f,
                                                        -3.0f, -3.0f))),
                        _mm256_mul_ps(xi_286, _mm256_set_ps(-3.0f, -3.0f, -3.0f,
                                                            -3.0f, -3.0f, -3.0f,
                                                            -3.0f, -3.0f))),
                    _mm256_mul_ps(xi_288,
                                  _mm256_set_ps(3.0f, 3.0f, 3.0f, 3.0f, 3.0f,
                                                3.0f, 3.0f, 3.0f))),
                xi_114),
            _mm256_set_ps(omega_even, omega_even, omega_even, omega_even,
                          omega_even, omega_even, omega_even, omega_even));
        const __m256 xi_116 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_mul_ps(xi_268,
                                  _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
                                                2.0f, 2.0f, 2.0f)),
                    _mm256_mul_ps(xi_280,
                                  _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
                                                2.0f, 2.0f, 2.0f))),
                _mm256_mul_ps(xi_284, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f,
                                                    2.0f, 2.0f, 2.0f, 2.0f))),
            _mm256_mul_ps(xi_286, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
                                                2.0f, 2.0f, 2.0f)));
        const __m256 xi_117 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_mul_ps(xi_273, _mm256_set_ps(5.0f, 5.0f, 5.0f, 5.0f,
                                                    5.0f, 5.0f, 5.0f, 5.0f)),
                _mm256_mul_ps(xi_282, _mm256_set_ps(5.0f, 5.0f, 5.0f, 5.0f,
                                                    5.0f, 5.0f, 5.0f, 5.0f))),
            xi_116);
        const __m256 xi_118 = _mm256_mul_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_add_ps(
                                        _mm256_mul_ps(
                                            xi_277,
                                            _mm256_set_ps(-2.0f, -2.0f, -2.0f,
                                                          -2.0f, -2.0f, -2.0f,
                                                          -2.0f, -2.0f)),
                                        _mm256_mul_ps(
                                            xi_281,
                                            _mm256_set_ps(-5.0f, -5.0f, -5.0f,
                                                          -5.0f, -5.0f, -5.0f,
                                                          -5.0f, -5.0f))),
                                    _mm256_mul_ps(xi_283,
                                                  _mm256_set_ps(-5.0f, -5.0f,
                                                                -5.0f, -5.0f,
                                                                -5.0f, -5.0f,
                                                                -5.0f, -5.0f))),
                                _mm256_mul_ps(xi_287,
                                              _mm256_set_ps(-5.0f, -5.0f, -5.0f,
                                                            -5.0f, -5.0f, -5.0f,
                                                            -5.0f, -5.0f))),
                            _mm256_mul_ps(xi_288,
                                          _mm256_set_ps(-2.0f, -2.0f, -2.0f,
                                                        -2.0f, -2.0f, -2.0f,
                                                        -2.0f, -2.0f))),
                        _mm256_mul_ps(xi_289, _mm256_set_ps(-5.0f, -5.0f, -5.0f,
                                                            -5.0f, -5.0f, -5.0f,
                                                            -5.0f, -5.0f))),
                    xi_114),
                xi_117),
            _mm256_set_ps(omega_even, omega_even, omega_even, omega_even,
                          omega_even, omega_even, omega_even, omega_even));
        const __m256 xi_121 =
            _mm256_mul_ps(xi_280, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_122 = _mm256_add_ps(xi_121, xi_19);
        const __m256 xi_123 =
            _mm256_mul_ps(xi_275, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_126 =
            _mm256_mul_ps(xi_283, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_127 = _mm256_add_ps(
            _mm256_add_ps(_mm256_add_ps(xi_12, xi_126), xi_16), xi_22);
        const __m256 xi_129 =
            _mm256_mul_ps(xi_289, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
                                                2.0f, 2.0f, 2.0f));
        const __m256 xi_130 =
            _mm256_mul_ps(xi_283, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
                                                2.0f, 2.0f, 2.0f));
        const __m256 xi_131 = _mm256_add_ps(
            _mm256_mul_ps(xi_281, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
                                                2.0f, 2.0f, 2.0f)),
            _mm256_mul_ps(xi_287, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
                                                2.0f, 2.0f, 2.0f)));
        const __m256 xi_132 = _mm256_mul_ps(
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
                                                                xi_270,
                                                                _mm256_set_ps(
                                                                    -4.0f,
                                                                    -4.0f,
                                                                    -4.0f,
                                                                    -4.0f,
                                                                    -4.0f,
                                                                    -4.0f,
                                                                    -4.0f,
                                                                    -4.0f)),
                                                            _mm256_mul_ps(
                                                                xi_272,
                                                                _mm256_set_ps(
                                                                    -4.0f,
                                                                    -4.0f,
                                                                    -4.0f,
                                                                    -4.0f,
                                                                    -4.0f,
                                                                    -4.0f,
                                                                    -4.0f,
                                                                    -4.0f))),
                                                        _mm256_mul_ps(
                                                            xi_275,
                                                            _mm256_set_ps(
                                                                -7.0f, -7.0f,
                                                                -7.0f, -7.0f,
                                                                -7.0f, -7.0f,
                                                                -7.0f, -7.0f))),
                                                    _mm256_mul_ps(
                                                        xi_277,
                                                        _mm256_set_ps(
                                                            5.0f, 5.0f, 5.0f,
                                                            5.0f, 5.0f, 5.0f,
                                                            5.0f, 5.0f))),
                                                _mm256_mul_ps(
                                                    xi_278,
                                                    _mm256_set_ps(
                                                        -7.0f, -7.0f, -7.0f,
                                                        -7.0f, -7.0f, -7.0f,
                                                        -7.0f, -7.0f))),
                                            _mm256_mul_ps(
                                                xi_279, _mm256_set_ps(
                                                            -7.0f, -7.0f, -7.0f,
                                                            -7.0f, -7.0f, -7.0f,
                                                            -7.0f, -7.0f))),
                                        _mm256_mul_ps(
                                            xi_285,
                                            _mm256_set_ps(-7.0f, -7.0f, -7.0f,
                                                          -7.0f, -7.0f, -7.0f,
                                                          -7.0f, -7.0f))),
                                    _mm256_mul_ps(
                                        xi_288,
                                        _mm256_set_ps(5.0f, 5.0f, 5.0f, 5.0f,
                                                      5.0f, 5.0f, 5.0f, 5.0f))),
                                xi_113),
                            xi_117),
                        xi_129),
                    xi_130),
                xi_131),
            _mm256_set_ps(omega_even, omega_even, omega_even, omega_even,
                          omega_even, omega_even, omega_even, omega_even));
        const __m256 xi_133 = _mm256_add_ps(xi_121, xi_284);
        const __m256 xi_134 = _mm256_add_ps(
            _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_133, xi_15), xi_23),
                          xi_286),
            xi_288);
        const __m256 xi_136 = _mm256_mul_ps(
            xi_134, _mm256_set_ps(xi_135, xi_135, xi_135, xi_135, xi_135,
                                  xi_135, xi_135, xi_135));
        const __m256 xi_138 = _mm256_add_ps(xi_123, xi_285);
        const __m256d xi_142 = _mm256_add_ps(
            random_5_1, _mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f,
                                      -0.5f, -0.5f));
        const __m256 xi_147 =
            _mm256_mul_ps(xi_278, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
                                                2.0f, 2.0f, 2.0f));
        const __m256 xi_148 =
            _mm256_mul_ps(xi_285, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
                                                2.0f, 2.0f, 2.0f));
        const __m256 xi_149 = _mm256_add_ps(
            _mm256_mul_ps(xi_275, _mm256_set_ps(-2.0f, -2.0f, -2.0f, -2.0f,
                                                -2.0f, -2.0f, -2.0f, -2.0f)),
            _mm256_mul_ps(xi_279, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f,
                                                2.0f, 2.0f, 2.0f)));
        const __m256 xi_150 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_mul_ps(
                                xi_147, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                            xi_148),
                        xi_149),
                    xi_15),
                xi_2),
            xi_20);
        const __m256 xi_152 = _mm256_mul_ps(
            xi_150, _mm256_set_ps(xi_151, xi_151, xi_151, xi_151, xi_151,
                                  xi_151, xi_151, xi_151));
        const __m256d xi_153 = _mm256_add_ps(
            random_3_0, _mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f,
                                      -0.5f, -0.5f));
        const __m256d xi_158 = _mm256_add_ps(
            random_0_1, _mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f,
                                      -0.5f, -0.5f));
        const __m256 xi_162 = _mm256_add_ps(xi_287, xi_289);
        const __m256 xi_176 = _mm256_add_ps(xi_126, xi_289);
        const __m256 xi_177 = _mm256_add_ps(
            _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_13, xi_176), xi_21),
                          xi_282),
            xi_287);
        const __m256 xi_178 = _mm256_mul_ps(
            xi_177, _mm256_set_ps(xi_135, xi_135, xi_135, xi_135, xi_135,
                                  xi_135, xi_135, xi_135));
        const __m256d xi_179 = _mm256_add_ps(
            random_4_1, _mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f,
                                      -0.5f, -0.5f));
        const __m256 xi_181 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_mul_ps(xi_148,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        xi_1),
                    xi_14),
                xi_147),
            xi_149);
        const __m256 xi_182 = _mm256_mul_ps(
            xi_181, _mm256_set_ps(xi_151, xi_151, xi_151, xi_151, xi_151,
                                  xi_151, xi_151, xi_151));
        const __m256d xi_183 = _mm256_add_ps(
            random_4_0, _mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f,
                                      -0.5f, -0.5f));
        const __m256 xi_188 = _mm256_add_ps(xi_268, xi_286);
        const __m256 xi_189 = _mm256_add_ps(
            _mm256_add_ps(_mm256_add_ps(xi_122, xi_188), xi_24), xi_272);
        const __m256 xi_190 = _mm256_mul_ps(
            xi_189, _mm256_set_ps(xi_135, xi_135, xi_135, xi_135, xi_135,
                                  xi_135, xi_135, xi_135));
        const __m256d xi_193 = _mm256_add_ps(
            random_5_0, _mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f,
                                      -0.5f, -0.5f));
        const __m256 xi_195 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_mul_ps(xi_129,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_130,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0))),
                    xi_131),
                xi_25),
            xi_3);
        const __m256 xi_196 = _mm256_mul_ps(
            xi_195, _mm256_set_ps(xi_151, xi_151, xi_151, xi_151, xi_151,
                                  xi_151, xi_151, xi_151));
        const __m256d xi_197 = _mm256_add_ps(
            random_3_1, _mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f,
                                      -0.5f, -0.5f));
        const __m256 xi_204 = _mm256_mul_ps(
            xi_132, _mm256_set_ps(0.0138888888888889f, 0.0138888888888889f,
                                  0.0138888888888889f, 0.0138888888888889f,
                                  0.0138888888888889f, 0.0138888888888889f,
                                  0.0138888888888889f, 0.0138888888888889f));
        const __m256 xi_225 = _mm256_mul_ps(
            xi_118,
            _mm256_set_ps(-0.00714285714285714f, -0.00714285714285714f,
                          -0.00714285714285714f, -0.00714285714285714f,
                          -0.00714285714285714f, -0.00714285714285714f,
                          -0.00714285714285714f, -0.00714285714285714f));
        const __m256 xi_227 = _mm256_mul_ps(
            xi_115, _mm256_set_ps(0.025f, 0.025f, 0.025f, 0.025f, 0.025f,
                                  0.025f, 0.025f, 0.025f));
        const __m256 xi_232 = _mm256_mul_ps(
            xi_195, _mm256_set_ps(xi_231, xi_231, xi_231, xi_231, xi_231,
                                  xi_231, xi_231, xi_231));
        const __m256 xi_234 = _mm256_mul_ps(
            xi_189, _mm256_set_ps(xi_233, xi_233, xi_233, xi_233, xi_233,
                                  xi_233, xi_233, xi_233));
        const __m256 xi_243 = _mm256_mul_ps(
            xi_150, _mm256_set_ps(xi_231, xi_231, xi_231, xi_231, xi_231,
                                  xi_231, xi_231, xi_231));
        const __m256 xi_244 = _mm256_mul_ps(
            xi_134, _mm256_set_ps(xi_233, xi_233, xi_233, xi_233, xi_233,
                                  xi_233, xi_233, xi_233));
        const __m256 xi_252 = _mm256_mul_ps(
            xi_118, _mm256_set_ps(0.0178571428571429f, 0.0178571428571429f,
                                  0.0178571428571429f, 0.0178571428571429f,
                                  0.0178571428571429f, 0.0178571428571429f,
                                  0.0178571428571429f, 0.0178571428571429f));
        const __m256 xi_258 = _mm256_mul_ps(
            xi_177, _mm256_set_ps(xi_233, xi_233, xi_233, xi_233, xi_233,
                                  xi_233, xi_233, xi_233));
        const __m256 xi_259 = _mm256_mul_ps(
            xi_181, _mm256_set_ps(xi_231, xi_231, xi_231, xi_231, xi_231,
                                  xi_231, xi_231, xi_231));
        const __m256 xi_32 =
            _mm256_mul_ps(xi_31, _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0,
                                               rr_0, rr_0, rr_0));
        const __m256 xi_45 =
            _mm256_mul_ps(xi_44, _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0,
                                               rr_0, rr_0, rr_0));
        const __m256 xi_52 =
            _mm256_mul_ps(xi_51, _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0,
                                               rr_0, rr_0, rr_0));
        const __m256 xi_55 =
            _mm256_mul_ps(xi_276, _mm256_set_ps(xi_54, xi_54, xi_54, xi_54,
                                                xi_54, xi_54, xi_54, xi_54));
        const __m256 xi_60 =
            _mm256_mul_ps(xi_271, _mm256_set_ps(xi_54, xi_54, xi_54, xi_54,
                                                xi_54, xi_54, xi_54, xi_54));
        const __m256 xi_82 =
            _mm256_mul_ps(xi_269, _mm256_set_ps(xi_54, xi_54, xi_54, xi_54,
                                                xi_54, xi_54, xi_54, xi_54));
        const __m256 vel0Term =
            _mm256_add_ps(_mm256_add_ps(xi_1, xi_275), xi_285);
        const __m256 vel1Term = _mm256_add_ps(xi_2, xi_278);
        const __m256 vel2Term = _mm256_add_ps(xi_289, xi_3);
        const __m256 rho = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_add_ps(vel0Term, vel1Term),
                                      vel2Term),
                        xi_274),
                    xi_4),
                xi_5),
            xi_6);
        const __m256 xi_7 = _mm256_div_ps(
            _mm256_set_ps(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0), rho);
        const __m256 xi_8 =
            _mm256_mul_ps(xi_7, _mm256_set_ps(0.5f, 0.5f, 0.5f, 0.5f, 0.5f,
                                              0.5f, 0.5f, 0.5f));
        const __m256 xi_106 =
            _mm256_mul_ps(rho, _mm256_set_ps(kT, kT, kT, kT, kT, kT, kT, kT));
        const __m256 xi_107 = _mm256_sqrt_ps(_mm256_mul_ps(
            xi_106,
            _mm256_set_ps(-((-omega_even + 1.0f) * (-omega_even + 1.0f)) + 1.0f,
                          -((-omega_even + 1.0f) * (-omega_even + 1.0f)) + 1.0f,
                          -((-omega_even + 1.0f) * (-omega_even + 1.0f)) + 1.0f,
                          -((-omega_even + 1.0f) * (-omega_even + 1.0f)) + 1.0f,
                          -((-omega_even + 1.0f) * (-omega_even + 1.0f)) + 1.0f,
                          -((-omega_even + 1.0f) * (-omega_even + 1.0f)) + 1.0f,
                          -((-omega_even + 1.0f) * (-omega_even + 1.0f)) + 1.0f,
                          -((-omega_even + 1.0f) * (-omega_even + 1.0f)) +
                              1.0f)));
        const __m256 xi_108 = _mm256_mul_ps(
            _mm256_mul_ps(
                _mm256_add_ps(random_6_0,
                              _mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f,
                                            -0.5f, -0.5f, -0.5f)),
                _mm256_set_ps(3.7416573867739413f, 3.7416573867739413f,
                              3.7416573867739413f, 3.7416573867739413f,
                              3.7416573867739413f, 3.7416573867739413f,
                              3.7416573867739413f, 3.7416573867739413f)),
            _mm256_set_ps(xi_107, xi_107, xi_107, xi_107, xi_107, xi_107,
                          xi_107, xi_107));
        const __m256 xi_109 = _mm256_mul_ps(
            _mm256_mul_ps(
                _mm256_add_ps(random_7_0,
                              _mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f,
                                            -0.5f, -0.5f, -0.5f)),
                _mm256_set_ps(5.4772255750516612f, 5.4772255750516612f,
                              5.4772255750516612f, 5.4772255750516612f,
                              5.4772255750516612f, 5.4772255750516612f,
                              5.4772255750516612f, 5.4772255750516612f)),
            _mm256_set_ps(xi_107, xi_107, xi_107, xi_107, xi_107, xi_107,
                          xi_107, xi_107));
        const __m256 xi_111 = _mm256_mul_ps(
            _mm256_mul_ps(
                _mm256_add_ps(random_2_1,
                              _mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f,
                                            -0.5f, -0.5f, -0.5f)),
                _mm256_set_ps(
                    _mm256_sqrt_ps(_mm256_mul_ps(
                        xi_106,
                        _mm256_set_ps(
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f))),
                    _mm256_sqrt_ps(_mm256_mul_ps(
                        xi_106,
                        _mm256_set_ps(
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f))),
                    _mm256_sqrt_ps(_mm256_mul_ps(
                        xi_106,
                        _mm256_set_ps(
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f))),
                    _mm256_sqrt_ps(_mm256_mul_ps(
                        xi_106,
                        _mm256_set_ps(
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f))),
                    _mm256_sqrt_ps(_mm256_mul_ps(
                        xi_106,
                        _mm256_set_ps(
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f))),
                    _mm256_sqrt_ps(_mm256_mul_ps(
                        xi_106,
                        _mm256_set_ps(
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f))),
                    _mm256_sqrt_ps(_mm256_mul_ps(
                        xi_106,
                        _mm256_set_ps(
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f))),
                    _mm256_sqrt_ps(_mm256_mul_ps(
                        xi_106,
                        _mm256_set_ps(
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f,
                            -((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) +
                                1.0f))))),
            _mm256_set_ps(xi_110, xi_110, xi_110, xi_110, xi_110, xi_110,
                          xi_110, xi_110));
        const __m256 xi_112 = _mm256_mul_ps(
            _mm256_mul_ps(
                _mm256_add_ps(random_6_1,
                              _mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f,
                                            -0.5f, -0.5f, -0.5f)),
                _mm256_set_ps(8.3666002653407556f, 8.3666002653407556f,
                              8.3666002653407556f, 8.3666002653407556f,
                              8.3666002653407556f, 8.3666002653407556f,
                              8.3666002653407556f, 8.3666002653407556f)),
            _mm256_set_ps(xi_107, xi_107, xi_107, xi_107, xi_107, xi_107,
                          xi_107, xi_107));
        const __m256 xi_143 = _mm256_sqrt_ps(_mm256_mul_ps(
            xi_106, _mm256_set_ps(
                        -((-omega_odd + 1.0f) * (-omega_odd + 1.0f)) + 1.0f,
                        -((-omega_odd + 1.0f) * (-omega_odd + 1.0f)) + 1.0f,
                        -((-omega_odd + 1.0f) * (-omega_odd + 1.0f)) + 1.0f,
                        -((-omega_odd + 1.0f) * (-omega_odd + 1.0f)) + 1.0f,
                        -((-omega_odd + 1.0f) * (-omega_odd + 1.0f)) + 1.0f,
                        -((-omega_odd + 1.0f) * (-omega_odd + 1.0f)) + 1.0f,
                        -((-omega_odd + 1.0f) * (-omega_odd + 1.0f)) + 1.0f,
                        -((-omega_odd + 1.0f) * (-omega_odd + 1.0f)) + 1.0f)));
        const __m256 xi_144 = _mm256_mul_ps(
            xi_143, _mm256_set_ps(1.4142135623730951f, 1.4142135623730951f,
                                  1.4142135623730951f, 1.4142135623730951f,
                                  1.4142135623730951f, 1.4142135623730951f,
                                  1.4142135623730951f, 1.4142135623730951f));
        const __m256 xi_145 =
            _mm256_mul_ps(xi_144, _mm256_set_ps(0.5f, 0.5f, 0.5f, 0.5f, 0.5f,
                                                0.5f, 0.5f, 0.5f));
        const __m256 xi_146 = _mm256_mul_ps(
            xi_142, _mm256_set_ps(xi_145, xi_145, xi_145, xi_145, xi_145,
                                  xi_145, xi_145, xi_145));
        const __m256 xi_154 = _mm256_mul_ps(
            xi_143, _mm256_set_ps(xi_110, xi_110, xi_110, xi_110, xi_110,
                                  xi_110, xi_110, xi_110));
        const __m256 xi_155 = _mm256_mul_ps(
            xi_154, _mm256_set_ps(0.166666666666667f, 0.166666666666667f,
                                  0.166666666666667f, 0.166666666666667f,
                                  0.166666666666667f, 0.166666666666667f,
                                  0.166666666666667f, 0.166666666666667f));
        const __m256 xi_156 = _mm256_mul_ps(
            xi_153, _mm256_set_ps(xi_155, xi_155, xi_155, xi_155, xi_155,
                                  xi_155, xi_155, xi_155));
        const __m256 xi_157 = _mm256_add_ps(
            _mm256_mul_ps(xi_152, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)),
            _mm256_mul_ps(xi_156, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)));
        const __m256 xi_159 = _mm256_sqrt_ps(_mm256_mul_ps(
            xi_106,
            _mm256_set_ps(
                -((-omega_shear + 1.0f) * (-omega_shear + 1.0f)) + 1.0f,
                -((-omega_shear + 1.0f) * (-omega_shear + 1.0f)) + 1.0f,
                -((-omega_shear + 1.0f) * (-omega_shear + 1.0f)) + 1.0f,
                -((-omega_shear + 1.0f) * (-omega_shear + 1.0f)) + 1.0f,
                -((-omega_shear + 1.0f) * (-omega_shear + 1.0f)) + 1.0f,
                -((-omega_shear + 1.0f) * (-omega_shear + 1.0f)) + 1.0f,
                -((-omega_shear + 1.0f) * (-omega_shear + 1.0f)) + 1.0f,
                -((-omega_shear + 1.0f) * (-omega_shear + 1.0f)) + 1.0f)));
        const __m256 xi_160 =
            _mm256_mul_ps(xi_159, _mm256_set_ps(0.5f, 0.5f, 0.5f, 0.5f, 0.5f,
                                                0.5f, 0.5f, 0.5f));
        const __m256 xi_161 = _mm256_mul_ps(
            xi_158, _mm256_set_ps(xi_160, xi_160, xi_160, xi_160, xi_160,
                                  xi_160, xi_160, xi_160));
        const __m256 xi_166 = _mm256_add_ps(
            _mm256_mul_ps(
                xi_108,
                _mm256_set_ps(-0.119047619047619f, -0.119047619047619f,
                              -0.119047619047619f, -0.119047619047619f,
                              -0.119047619047619f, -0.119047619047619f,
                              -0.119047619047619f, -0.119047619047619f)),
            _mm256_mul_ps(
                xi_132,
                _mm256_set_ps(-0.0198412698412698f, -0.0198412698412698f,
                              -0.0198412698412698f, -0.0198412698412698f,
                              -0.0198412698412698f, -0.0198412698412698f,
                              -0.0198412698412698f, -0.0198412698412698f)));
        const __m256 xi_168 = _mm256_mul_ps(
            _mm256_mul_ps(
                _mm256_add_ps(random_0_0,
                              _mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f,
                                            -0.5f, -0.5f, -0.5f)),
                _mm256_set_ps(1.7320508075688772f, 1.7320508075688772f,
                              1.7320508075688772f, 1.7320508075688772f,
                              1.7320508075688772f, 1.7320508075688772f,
                              1.7320508075688772f, 1.7320508075688772f)),
            _mm256_set_ps(xi_159, xi_159, xi_159, xi_159, xi_159, xi_159,
                          xi_159, xi_159));
        const __m256 xi_172 = _mm256_add_ps(xi_152, xi_156);
        const __m256 xi_180 = _mm256_mul_ps(
            xi_179, _mm256_set_ps(xi_145, xi_145, xi_145, xi_145, xi_145,
                                  xi_145, xi_145, xi_145));
        const __m256 xi_184 = _mm256_mul_ps(
            xi_183, _mm256_set_ps(xi_155, xi_155, xi_155, xi_155, xi_155,
                                  xi_155, xi_155, xi_155));
        const __m256 xi_185 = _mm256_add_ps(xi_182, xi_184);
        const __m256 xi_187 = _mm256_add_ps(
            _mm256_mul_ps(xi_182, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)),
            _mm256_mul_ps(xi_184, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)));
        const __m256 xi_194 = _mm256_mul_ps(
            xi_193, _mm256_set_ps(xi_145, xi_145, xi_145, xi_145, xi_145,
                                  xi_145, xi_145, xi_145));
        const __m256 xi_198 = _mm256_mul_ps(
            xi_197, _mm256_set_ps(xi_155, xi_155, xi_155, xi_155, xi_155,
                                  xi_155, xi_155, xi_155));
        const __m256 xi_199 = _mm256_add_ps(
            _mm256_mul_ps(xi_196, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)),
            _mm256_mul_ps(xi_198, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)));
        const __m256 xi_201 = _mm256_add_ps(xi_196, xi_198);
        const __m256 xi_202 = _mm256_mul_ps(
            _mm256_mul_ps(xi_158, _mm256_set_ps(0.25f, 0.25f, 0.25f, 0.25f,
                                                0.25f, 0.25f, 0.25f, 0.25f)),
            _mm256_set_ps(xi_159, xi_159, xi_159, xi_159, xi_159, xi_159,
                          xi_159, xi_159));
        const __m256 xi_205 = _mm256_mul_ps(
            xi_108, _mm256_set_ps(0.0833333333333333f, 0.0833333333333333f,
                                  0.0833333333333333f, 0.0833333333333333f,
                                  0.0833333333333333f, 0.0833333333333333f,
                                  0.0833333333333333f, 0.0833333333333333f));
        const __m256 xi_215 = _mm256_mul_ps(
            _mm256_add_ps(random_1_0,
                          _mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f,
                                        -0.5f, -0.5f, -0.5f)),
            _mm256_set_ps(xi_160, xi_160, xi_160, xi_160, xi_160, xi_160,
                          xi_160, xi_160));
        const __m256 xi_224 = _mm256_mul_ps(
            _mm256_add_ps(random_2_0,
                          _mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f,
                                        -0.5f, -0.5f, -0.5f)),
            _mm256_set_ps(xi_160, xi_160, xi_160, xi_160, xi_160, xi_160,
                          xi_160, xi_160));
        const __m256 xi_228 = _mm256_mul_ps(
            xi_112, _mm256_set_ps(-0.0142857142857143f, -0.0142857142857143f,
                                  -0.0142857142857143f, -0.0142857142857143f,
                                  -0.0142857142857143f, -0.0142857142857143f,
                                  -0.0142857142857143f, -0.0142857142857143f));
        const __m256 xi_229 =
            _mm256_mul_ps(xi_109, _mm256_set_ps(0.05f, 0.05f, 0.05f, 0.05f,
                                                0.05f, 0.05f, 0.05f, 0.05f));
        const __m256 xi_235 = _mm256_mul_ps(
            xi_154, _mm256_set_ps(0.0833333333333333f, 0.0833333333333333f,
                                  0.0833333333333333f, 0.0833333333333333f,
                                  0.0833333333333333f, 0.0833333333333333f,
                                  0.0833333333333333f, 0.0833333333333333f));
        const __m256 xi_236 = _mm256_mul_ps(
            xi_197, _mm256_set_ps(xi_235, xi_235, xi_235, xi_235, xi_235,
                                  xi_235, xi_235, xi_235));
        const __m256 xi_237 =
            _mm256_mul_ps(xi_144, _mm256_set_ps(0.25f, 0.25f, 0.25f, 0.25f,
                                                0.25f, 0.25f, 0.25f, 0.25f));
        const __m256 xi_238 = _mm256_mul_ps(
            xi_193, _mm256_set_ps(xi_237, xi_237, xi_237, xi_237, xi_237,
                                  xi_237, xi_237, xi_237));
        const __m256 xi_240 = _mm256_add_ps(
            _mm256_mul_ps(
                xi_108,
                _mm256_set_ps(-0.0238095238095238f, -0.0238095238095238f,
                              -0.0238095238095238f, -0.0238095238095238f,
                              -0.0238095238095238f, -0.0238095238095238f,
                              -0.0238095238095238f, -0.0238095238095238f)),
            _mm256_mul_ps(
                xi_132,
                _mm256_set_ps(-0.00396825396825397f, -0.00396825396825397f,
                              -0.00396825396825397f, -0.00396825396825397f,
                              -0.00396825396825397f, -0.00396825396825397f,
                              -0.00396825396825397f, -0.00396825396825397f)));
        const __m256 xi_245 = _mm256_mul_ps(
            xi_153, _mm256_set_ps(xi_235, xi_235, xi_235, xi_235, xi_235,
                                  xi_235, xi_235, xi_235));
        const __m256 xi_246 = _mm256_mul_ps(
            xi_142, _mm256_set_ps(xi_237, xi_237, xi_237, xi_237, xi_237,
                                  xi_237, xi_237, xi_237));
        const __m256 xi_250 =
            _mm256_mul_ps(xi_202, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_253 = _mm256_mul_ps(
            xi_112, _mm256_set_ps(0.0357142857142857f, 0.0357142857142857f,
                                  0.0357142857142857f, 0.0357142857142857f,
                                  0.0357142857142857f, 0.0357142857142857f,
                                  0.0357142857142857f, 0.0357142857142857f));
        const __m256 xi_255 = _mm256_mul_ps(
            _mm256_add_ps(random_1_1,
                          _mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f,
                                        -0.5f, -0.5f, -0.5f)),
            _mm256_set_ps(xi_160, xi_160, xi_160, xi_160, xi_160, xi_160,
                          xi_160, xi_160));
        const __m256 xi_260 = _mm256_mul_ps(
            xi_179, _mm256_set_ps(xi_237, xi_237, xi_237, xi_237, xi_237,
                                  xi_237, xi_237, xi_237));
        const __m256 xi_261 = _mm256_mul_ps(
            xi_183, _mm256_set_ps(xi_235, xi_235, xi_235, xi_235, xi_235,
                                  xi_235, xi_235, xi_235));
        const __m256 u_0 = _mm256_add_ps(
            _mm256_mul_ps(xi_7,
                          _mm256_add_ps(_mm256_add_ps(vel0Term, xi_10), xi_14)),
            _mm256_mul_ps(xi_276, xi_8));
        const __m256 xi_26 = _mm256_mul_ps(u_0, xi_276);
        const __m256 xi_38 = _mm256_mul_ps(
            xi_26, _mm256_set_ps(0.166666666666667f, 0.166666666666667f,
                                 0.166666666666667f, 0.166666666666667f,
                                 0.166666666666667f, 0.166666666666667f,
                                 0.166666666666667f, 0.166666666666667f));
        const __m256 xi_39 = _mm256_mul_ps(
            xi_26, _mm256_set_ps(0.0833333333333333f, 0.0833333333333333f,
                                 0.0833333333333333f, 0.0833333333333333f,
                                 0.0833333333333333f, 0.0833333333333333f,
                                 0.0833333333333333f, 0.0833333333333333f));
        const __m256 xi_40 = _mm256_mul_ps(
            xi_39,
            _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear,
                          omega_shear, omega_shear, omega_shear, omega_shear));
        const __m256 xi_41 = _mm256_add_ps(
            _mm256_mul_ps(xi_38, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                               -1.0, -1.0, -1.0)),
            xi_40);
        const __m256 xi_57 = _mm256_add_ps(
            _mm256_mul_ps(
                _mm256_mul_ps(xi_26, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_set_ps(xi_56, xi_56, xi_56, xi_56, xi_56, xi_56, xi_56,
                              xi_56)),
            xi_38);
        const __m256 xi_58 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_mul_ps(xi_44, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                xi_55),
            xi_57);
        const __m256 xi_62 = _mm256_mul_ps(
            _mm256_mul_ps(xi_26, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                               -1.0, -1.0, -1.0)),
            _mm256_set_ps(xi_61, xi_61, xi_61, xi_61, xi_61, xi_61, xi_61,
                          xi_61));
        const __m256 xi_69 = _mm256_mul_ps(u_0, xi_68);
        const __m256 xi_74 = _mm256_mul_ps(u_0, xi_73);
        const __m256 xi_78 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_mul_ps(xi_55, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                xi_44),
            xi_57);
        const __m256 xi_85 =
            _mm256_mul_ps(xi_39, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                               -1.0, -1.0, -1.0));
        const __m256 xi_96 = _mm256_mul_ps(u_0, xi_269);
        const __m256 xi_97 =
            _mm256_mul_ps(xi_96, _mm256_set_ps(0.25f, 0.25f, 0.25f, 0.25f,
                                               0.25f, 0.25f, 0.25f, 0.25f));
        const __m256 xi_100 =
            _mm256_mul_ps(xi_96, _mm256_set_ps(xi_72, xi_72, xi_72, xi_72,
                                               xi_72, xi_72, xi_72, xi_72));
        const __m256 xi_119 = _mm256_mul_ps(rho, (_mm256_mul_ps(u_0, u_0)));
        const __m256 xi_173 = _mm256_mul_ps(rho, u_0);
        const __m256 xi_174 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_mul_ps(vel0Term,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        xi_162),
                    xi_173),
                xi_278),
            xi_4);
        const __m256 xi_175 = _mm256_mul_ps(
            xi_174, _mm256_set_ps(xi_140, xi_140, xi_140, xi_140, xi_140,
                                  xi_140, xi_140, xi_140));
        const __m256 xi_211 = _mm256_mul_ps(
            xi_174, _mm256_set_ps(xi_206, xi_206, xi_206, xi_206, xi_206,
                                  xi_206, xi_206, xi_206));
        const __m256 u_1 = _mm256_add_ps(
            _mm256_mul_ps(
                xi_7,
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_add_ps(vel1Term, xi_17), xi_20),
                        xi_275),
                    xi_9)),
            _mm256_mul_ps(xi_271, xi_8));
        const __m256 xi_27 = _mm256_mul_ps(u_1, xi_271);
        const __m256 xi_33 = _mm256_mul_ps(
            xi_27, _mm256_set_ps(0.166666666666667f, 0.166666666666667f,
                                 0.166666666666667f, 0.166666666666667f,
                                 0.166666666666667f, 0.166666666666667f,
                                 0.166666666666667f, 0.166666666666667f));
        const __m256 xi_46 = _mm256_mul_ps(
            xi_27, _mm256_set_ps(0.0833333333333333f, 0.0833333333333333f,
                                 0.0833333333333333f, 0.0833333333333333f,
                                 0.0833333333333333f, 0.0833333333333333f,
                                 0.0833333333333333f, 0.0833333333333333f));
        const __m256 xi_47 = _mm256_mul_ps(
            xi_46,
            _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear,
                          omega_shear, omega_shear, omega_shear, omega_shear));
        const __m256 xi_48 = _mm256_add_ps(
            _mm256_mul_ps(xi_33, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                               -1.0, -1.0, -1.0)),
            xi_47);
        const __m256 xi_63 = _mm256_mul_ps(
            _mm256_mul_ps(xi_27, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                               -1.0, -1.0, -1.0)),
            _mm256_set_ps(xi_61, xi_61, xi_61, xi_61, xi_61, xi_61, xi_61,
                          xi_61));
        const __m256 xi_70 =
            _mm256_mul_ps(u_1, _mm256_set_ps(0.25f, 0.25f, 0.25f, 0.25f, 0.25f,
                                             0.25f, 0.25f, 0.25f));
        const __m256 xi_71 = _mm256_mul_ps(xi_276, xi_70);
        const __m256 xi_75 =
            _mm256_mul_ps(u_1, _mm256_set_ps(xi_72, xi_72, xi_72, xi_72, xi_72,
                                             xi_72, xi_72, xi_72));
        const __m256 xi_76 = _mm256_mul_ps(xi_276, xi_75);
        const __m256 xi_77 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_mul_ps(xi_69, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                       -1.0, -1.0, -1.0, -1.0)),
                    _mm256_mul_ps(xi_71,
                                  _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0))),
                xi_74),
            xi_76);
        const __m256 xi_79 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_mul_ps(xi_74, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                       -1.0, -1.0, -1.0, -1.0)),
                    _mm256_mul_ps(xi_76,
                                  _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0))),
                xi_69),
            xi_71);
        const __m256 xi_87 = _mm256_mul_ps(xi_269, xi_70);
        const __m256 xi_89 = _mm256_mul_ps(xi_269, xi_75);
        const __m256 xi_94 =
            _mm256_mul_ps(xi_46, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                               -1.0, -1.0, -1.0));
        const __m256 xi_124 = _mm256_mul_ps(rho, (_mm256_mul_ps(u_1, u_1)));
        const __m256 xi_125 =
            _mm256_add_ps(_mm256_add_ps(xi_10, xi_123), xi_124);
        const __m256 xi_137 = _mm256_mul_ps(rho, u_1);
        const __m256 xi_139 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_mul_ps(vel1Term,
                                          _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                        -1.0, -1.0, -1.0,
                                                        -1.0)),
                            xi_137),
                        xi_138),
                    xi_279),
                xi_284),
            xi_5);
        const __m256 xi_141 = _mm256_mul_ps(
            xi_139, _mm256_set_ps(xi_140, xi_140, xi_140, xi_140, xi_140,
                                  xi_140, xi_140, xi_140));
        const __m256 xi_207 = _mm256_mul_ps(
            xi_139, _mm256_set_ps(xi_206, xi_206, xi_206, xi_206, xi_206,
                                  xi_206, xi_206, xi_206));
        const __m256 xi_217 = _mm256_mul_ps(
            _mm256_add_ps(
                _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(u_0, xi_137), xi_138),
                              xi_278),
                xi_9),
            _mm256_set_ps(xi_216, xi_216, xi_216, xi_216, xi_216, xi_216,
                          xi_216, xi_216));
        const __m256 xi_218 = _mm256_add_ps(
            _mm256_mul_ps(xi_215, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)),
            _mm256_mul_ps(xi_217, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)));
        const __m256 xi_219 = _mm256_add_ps(xi_215, xi_217);
        const __m256 u_2 = _mm256_add_ps(
            _mm256_mul_ps(
                xi_7, _mm256_add_ps(
                          _mm256_add_ps(_mm256_add_ps(vel2Term, xi_22), xi_25),
                          xi_283)),
            _mm256_mul_ps(xi_269, xi_8));
        const __m256 xi_28 = _mm256_mul_ps(u_2, xi_269);
        const __m256 xi_34 = _mm256_mul_ps(
            xi_28, _mm256_set_ps(0.166666666666667f, 0.166666666666667f,
                                 0.166666666666667f, 0.166666666666667f,
                                 0.166666666666667f, 0.166666666666667f,
                                 0.166666666666667f, 0.166666666666667f));
        const __m256 xi_35 = _mm256_mul_ps(
            xi_28, _mm256_set_ps(0.0833333333333333f, 0.0833333333333333f,
                                 0.0833333333333333f, 0.0833333333333333f,
                                 0.0833333333333333f, 0.0833333333333333f,
                                 0.0833333333333333f, 0.0833333333333333f));
        const __m256 xi_36 = _mm256_mul_ps(
            xi_35,
            _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear,
                          omega_shear, omega_shear, omega_shear, omega_shear));
        const __m256 xi_37 = _mm256_add_ps(
            _mm256_mul_ps(xi_34, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                               -1.0, -1.0, -1.0)),
            xi_36);
        const __m256 xi_42 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_mul_ps(
                        xi_27,
                        _mm256_set_ps(0.333333333333333f, 0.333333333333333f,
                                      0.333333333333333f, 0.333333333333333f,
                                      0.333333333333333f, 0.333333333333333f,
                                      0.333333333333333f, 0.333333333333333f)),
                    _mm256_mul_ps(
                        _mm256_mul_ps(xi_33,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                      omega_shear, omega_shear, omega_shear,
                                      omega_shear, omega_shear))),
                xi_37),
            xi_41);
        const __m256 xi_49 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_mul_ps(
                        xi_26,
                        _mm256_set_ps(0.333333333333333f, 0.333333333333333f,
                                      0.333333333333333f, 0.333333333333333f,
                                      0.333333333333333f, 0.333333333333333f,
                                      0.333333333333333f, 0.333333333333333f)),
                    _mm256_mul_ps(
                        _mm256_mul_ps(xi_38,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                      omega_shear, omega_shear, omega_shear,
                                      omega_shear, omega_shear))),
                xi_37),
            xi_48);
        const __m256 xi_53 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_mul_ps(
                        xi_28,
                        _mm256_set_ps(0.333333333333333f, 0.333333333333333f,
                                      0.333333333333333f, 0.333333333333333f,
                                      0.333333333333333f, 0.333333333333333f,
                                      0.333333333333333f, 0.333333333333333f)),
                    _mm256_mul_ps(
                        _mm256_mul_ps(xi_34,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_set_ps(omega_shear, omega_shear, omega_shear,
                                      omega_shear, omega_shear, omega_shear,
                                      omega_shear, omega_shear))),
                xi_41),
            xi_48);
        const __m256 xi_59 =
            _mm256_mul_ps(xi_35, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                               -1.0, -1.0, -1.0));
        const __m256 xi_64 = _mm256_mul_ps(
            _mm256_mul_ps(xi_28, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                               -1.0, -1.0, -1.0)),
            _mm256_set_ps(xi_61, xi_61, xi_61, xi_61, xi_61, xi_61, xi_61,
                          xi_61));
        const __m256 xi_65 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_mul_ps(
                            _mm256_mul_ps(xi_27, _mm256_set_ps(-1.0, -1.0, -1.0,
                                                               -1.0, -1.0, -1.0,
                                                               -1.0, -1.0)),
                            _mm256_set_ps(xi_56, xi_56, xi_56, xi_56, xi_56,
                                          xi_56, xi_56, xi_56)),
                        xi_33),
                    xi_62),
                xi_63),
            xi_64);
        const __m256 xi_66 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_mul_ps(xi_60, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                xi_31),
            xi_65);
        const __m256 xi_67 = _mm256_add_ps(_mm256_add_ps(xi_36, xi_59), xi_66);
        const __m256 xi_80 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_mul_ps(xi_31, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                xi_60),
            xi_65);
        const __m256 xi_81 = _mm256_add_ps(_mm256_add_ps(xi_36, xi_59), xi_80);
        const __m256 xi_83 = _mm256_add_ps(
            _mm256_mul_ps(
                _mm256_mul_ps(xi_28, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                _mm256_set_ps(xi_56, xi_56, xi_56, xi_56, xi_56, xi_56, xi_56,
                              xi_56)),
            xi_34);
        const __m256 xi_84 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_mul_ps(xi_82, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                xi_51),
            xi_83);
        const __m256 xi_86 = _mm256_add_ps(_mm256_add_ps(xi_40, xi_66), xi_85);
        const __m256 xi_88 = _mm256_mul_ps(u_2, xi_68);
        const __m256 xi_90 = _mm256_mul_ps(u_2, xi_73);
        const __m256 xi_91 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_mul_ps(xi_89, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                       -1.0, -1.0, -1.0, -1.0)),
                    _mm256_mul_ps(xi_90,
                                  _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0))),
                xi_87),
            xi_88);
        const __m256 xi_92 = _mm256_add_ps(_mm256_add_ps(xi_40, xi_80), xi_85);
        const __m256 xi_93 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_mul_ps(xi_87, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                       -1.0, -1.0, -1.0, -1.0)),
                    _mm256_mul_ps(xi_88,
                                  _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0))),
                xi_89),
            xi_90);
        const __m256 xi_95 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_47, xi_62), xi_63),
                              xi_64),
                xi_84),
            xi_94);
        const __m256 xi_98 = _mm256_mul_ps(u_2, xi_276);
        const __m256 xi_99 =
            _mm256_mul_ps(xi_98, _mm256_set_ps(0.25f, 0.25f, 0.25f, 0.25f,
                                               0.25f, 0.25f, 0.25f, 0.25f));
        const __m256 xi_101 =
            _mm256_mul_ps(xi_98, _mm256_set_ps(xi_72, xi_72, xi_72, xi_72,
                                               xi_72, xi_72, xi_72, xi_72));
        const __m256 xi_102 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_mul_ps(xi_97, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                       -1.0, -1.0, -1.0, -1.0)),
                    _mm256_mul_ps(xi_99,
                                  _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0))),
                xi_100),
            xi_101);
        const __m256 xi_103 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_mul_ps(xi_100,
                                  _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)),
                    _mm256_mul_ps(xi_101,
                                  _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0))),
                xi_97),
            xi_99);
        const __m256 xi_104 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_mul_ps(xi_51, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                xi_82),
            xi_83);
        const __m256 xi_105 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(xi_104, xi_47), xi_62), xi_63),
                xi_64),
            xi_94);
        const __m256 xi_120 = _mm256_mul_ps(rho, (_mm256_mul_ps(u_2, u_2)));
        const __m256 xi_128 = _mm256_mul_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(_mm256_add_ps(xi_119, xi_120),
                                              xi_122),
                                xi_125),
                            xi_127),
                        xi_18),
                    xi_23),
                xi_274),
            _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk, omega_bulk,
                          omega_bulk, omega_bulk, omega_bulk, omega_bulk));
        const __m256 xi_163 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_mul_ps(xi_120, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                xi_270),
            xi_272);
        const __m256 xi_164 = _mm256_mul_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_mul_ps(xi_288,
                                              _mm256_set_ps(-1.0, -1.0, -1.0,
                                                            -1.0, -1.0, -1.0,
                                                            -1.0, -1.0)),
                                xi_0),
                            xi_125),
                        xi_162),
                    xi_163),
                xi_17),
            _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear,
                          omega_shear, omega_shear, omega_shear, omega_shear));
        const __m256 xi_165 = _mm256_mul_ps(
            xi_164, _mm256_set_ps(0.125f, 0.125f, 0.125f, 0.125f, 0.125f,
                                  0.125f, 0.125f, 0.125f));
        const __m256 xi_167 = _mm256_mul_ps(
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
                                                        xi_119,
                                                        _mm256_set_ps(
                                                            2.0f, 2.0f, 2.0f,
                                                            2.0f, 2.0f, 2.0f,
                                                            2.0f, 2.0f)),
                                                    _mm256_mul_ps(
                                                        xi_124,
                                                        _mm256_set_ps(
                                                            -1.0, -1.0, -1.0,
                                                            -1.0, -1.0, -1.0,
                                                            -1.0, -1.0))),
                                                _mm256_mul_ps(
                                                    xi_273,
                                                    _mm256_set_ps(
                                                        -2.0f, -2.0f, -2.0f,
                                                        -2.0f, -2.0f, -2.0f,
                                                        -2.0f, -2.0f))),
                                            _mm256_mul_ps(
                                                xi_282, _mm256_set_ps(
                                                            -2.0f, -2.0f, -2.0f,
                                                            -2.0f, -2.0f, -2.0f,
                                                            -2.0f, -2.0f))),
                                        xi_10),
                                    xi_116),
                                xi_123),
                            xi_127),
                        xi_163),
                    xi_277),
                xi_288),
            _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear,
                          omega_shear, omega_shear, omega_shear, omega_shear));
        const __m256 xi_169 = _mm256_add_ps(
            _mm256_mul_ps(
                xi_167,
                _mm256_set_ps(-0.0416666666666667f, -0.0416666666666667f,
                              -0.0416666666666667f, -0.0416666666666667f,
                              -0.0416666666666667f, -0.0416666666666667f,
                              -0.0416666666666667f, -0.0416666666666667f)),
            _mm256_mul_ps(
                xi_168,
                _mm256_set_ps(-0.166666666666667f, -0.166666666666667f,
                              -0.166666666666667f, -0.166666666666667f,
                              -0.166666666666667f, -0.166666666666667f,
                              -0.166666666666667f, -0.166666666666667f)));
        const __m256 xi_170 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_mul_ps(xi_109,
                              _mm256_set_ps(-0.1f, -0.1f, -0.1f, -0.1f, -0.1f,
                                            -0.1f, -0.1f, -0.1f)),
                _mm256_mul_ps(xi_115,
                              _mm256_set_ps(-0.05f, -0.05f, -0.05f, -0.05f,
                                            -0.05f, -0.05f, -0.05f, -0.05f))),
            xi_169);
        const __m256 xi_171 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_mul_ps(
                                xi_112,
                                _mm256_set_ps(
                                    0.0285714285714286f, 0.0285714285714286f,
                                    0.0285714285714286f, 0.0285714285714286f,
                                    0.0285714285714286f, 0.0285714285714286f,
                                    0.0285714285714286f, 0.0285714285714286f)),
                            _mm256_mul_ps(
                                xi_118,
                                _mm256_set_ps(
                                    0.0142857142857143f, 0.0142857142857143f,
                                    0.0142857142857143f, 0.0142857142857143f,
                                    0.0142857142857143f, 0.0142857142857143f,
                                    0.0142857142857143f, 0.0142857142857143f))),
                        xi_161),
                    xi_165),
                xi_166),
            xi_170);
        const __m256 xi_186 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_mul_ps(
                            xi_112,
                            _mm256_set_ps(
                                -0.0714285714285714f, -0.0714285714285714f,
                                -0.0714285714285714f, -0.0714285714285714f,
                                -0.0714285714285714f, -0.0714285714285714f,
                                -0.0714285714285714f, -0.0714285714285714f)),
                        _mm256_mul_ps(
                            xi_118,
                            _mm256_set_ps(
                                -0.0357142857142857f, -0.0357142857142857f,
                                -0.0357142857142857f, -0.0357142857142857f,
                                -0.0357142857142857f, -0.0357142857142857f,
                                -0.0357142857142857f, -0.0357142857142857f))),
                    _mm256_mul_ps(
                        xi_167, _mm256_set_ps(
                                    0.0833333333333333f, 0.0833333333333333f,
                                    0.0833333333333333f, 0.0833333333333333f,
                                    0.0833333333333333f, 0.0833333333333333f,
                                    0.0833333333333333f, 0.0833333333333333f))),
                _mm256_mul_ps(
                    xi_168,
                    _mm256_set_ps(0.333333333333333f, 0.333333333333333f,
                                  0.333333333333333f, 0.333333333333333f,
                                  0.333333333333333f, 0.333333333333333f,
                                  0.333333333333333f, 0.333333333333333f))),
            xi_166);
        const __m256 xi_191 = _mm256_add_ps(
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
                            xi_121),
                        xi_126),
                    xi_188),
                xi_281),
            xi_6);
        const __m256 xi_192 = _mm256_mul_ps(
            xi_191, _mm256_set_ps(xi_140, xi_140, xi_140, xi_140, xi_140,
                                  xi_140, xi_140, xi_140));
        const __m256 xi_200 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_mul_ps(
                                    xi_108, _mm256_set_ps(0.0952380952380952f,
                                                          0.0952380952380952f,
                                                          0.0952380952380952f,
                                                          0.0952380952380952f,
                                                          0.0952380952380952f,
                                                          0.0952380952380952f,
                                                          0.0952380952380952f,
                                                          0.0952380952380952f)),
                                _mm256_mul_ps(
                                    xi_112,
                                    _mm256_set_ps(-0.0428571428571429f,
                                                  -0.0428571428571429f,
                                                  -0.0428571428571429f,
                                                  -0.0428571428571429f,
                                                  -0.0428571428571429f,
                                                  -0.0428571428571429f,
                                                  -0.0428571428571429f,
                                                  -0.0428571428571429f))),
                            _mm256_mul_ps(
                                xi_118,
                                _mm256_set_ps(
                                    -0.0214285714285714f, -0.0214285714285714f,
                                    -0.0214285714285714f, -0.0214285714285714f,
                                    -0.0214285714285714f, -0.0214285714285714f,
                                    -0.0214285714285714f,
                                    -0.0214285714285714f))),
                        _mm256_mul_ps(
                            xi_132,
                            _mm256_set_ps(
                                0.0158730158730159f, 0.0158730158730159f,
                                0.0158730158730159f, 0.0158730158730159f,
                                0.0158730158730159f, 0.0158730158730159f,
                                0.0158730158730159f, 0.0158730158730159f))),
                    _mm256_mul_ps(xi_161,
                                  _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0))),
                _mm256_mul_ps(xi_165, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0))),
            xi_170);
        const __m256 xi_203 = _mm256_mul_ps(
            xi_164, _mm256_set_ps(0.0625f, 0.0625f, 0.0625f, 0.0625f, 0.0625f,
                                  0.0625f, 0.0625f, 0.0625f));
        const __m256 xi_208 = _mm256_add_ps(
            _mm256_mul_ps(
                xi_111,
                _mm256_set_ps(0.0833333333333333f, 0.0833333333333333f,
                              0.0833333333333333f, 0.0833333333333333f,
                              0.0833333333333333f, 0.0833333333333333f,
                              0.0833333333333333f, 0.0833333333333333f)),
            _mm256_mul_ps(
                xi_128,
                _mm256_set_ps(0.0416666666666667f, 0.0416666666666667f,
                              0.0416666666666667f, 0.0416666666666667f,
                              0.0416666666666667f, 0.0416666666666667f,
                              0.0416666666666667f, 0.0416666666666667f)));
        const __m256 xi_209 = _mm256_add_ps(xi_207, xi_208);
        const __m256 xi_210 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(xi_172, xi_202), xi_203),
                    xi_204),
                xi_205),
            xi_209);
        const __m256 xi_212 = _mm256_add_ps(
            _mm256_mul_ps(
                xi_167,
                _mm256_set_ps(0.0208333333333333f, 0.0208333333333333f,
                              0.0208333333333333f, 0.0208333333333333f,
                              0.0208333333333333f, 0.0208333333333333f,
                              0.0208333333333333f, 0.0208333333333333f)),
            _mm256_mul_ps(
                xi_168,
                _mm256_set_ps(0.0833333333333333f, 0.0833333333333333f,
                              0.0833333333333333f, 0.0833333333333333f,
                              0.0833333333333333f, 0.0833333333333333f,
                              0.0833333333333333f, 0.0833333333333333f)));
        const __m256 xi_213 = _mm256_add_ps(
            _mm256_mul_ps(xi_211, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)),
            xi_212);
        const __m256 xi_214 = _mm256_add_ps(xi_187, xi_213);
        const __m256 xi_220 = _mm256_add_ps(xi_211, xi_212);
        const __m256 xi_221 = _mm256_add_ps(xi_185, xi_220);
        const __m256 xi_222 = _mm256_add_ps(
            _mm256_mul_ps(xi_207, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)),
            xi_208);
        const __m256 xi_223 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(xi_157, xi_202), xi_203),
                    xi_204),
                xi_205),
            xi_222);
        const __m256 xi_226 = _mm256_mul_ps(
            _mm256_add_ps(
                _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(u_2, xi_137), xi_133),
                              xi_18),
                xi_268),
            _mm256_set_ps(xi_216, xi_216, xi_216, xi_216, xi_216, xi_216,
                          xi_216, xi_216));
        const __m256 xi_230 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_add_ps(xi_169, xi_224), xi_225),
                        xi_226),
                    xi_227),
                xi_228),
            xi_229);
        const __m256 xi_239 = _mm256_mul_ps(
            xi_191, _mm256_set_ps(xi_206, xi_206, xi_206, xi_206, xi_206,
                                  xi_206, xi_206, xi_206));
        const __m256 xi_241 = _mm256_add_ps(xi_239, xi_240);
        const __m256 xi_242 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_mul_ps(xi_232,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_236,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0))),
                    xi_234),
                xi_238),
            xi_241);
        const __m256 xi_247 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_mul_ps(xi_243,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_245,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0))),
                    xi_209),
                xi_244),
            xi_246);
        const __m256 xi_248 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_mul_ps(xi_244,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_246,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0))),
                    xi_222),
                xi_243),
            xi_245);
        const __m256 xi_249 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_mul_ps(xi_224,
                                              _mm256_set_ps(-1.0, -1.0, -1.0,
                                                            -1.0, -1.0, -1.0,
                                                            -1.0, -1.0)),
                                _mm256_mul_ps(xi_226,
                                              _mm256_set_ps(-1.0, -1.0, -1.0,
                                                            -1.0, -1.0, -1.0,
                                                            -1.0, -1.0))),
                            xi_169),
                        xi_225),
                    xi_227),
                xi_228),
            xi_229);
        const __m256 xi_251 =
            _mm256_mul_ps(xi_203, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0));
        const __m256 xi_254 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_add_ps(xi_201, xi_208), xi_241),
                        xi_250),
                    xi_251),
                xi_252),
            xi_253);
        const __m256 xi_256 = _mm256_mul_ps(
            _mm256_add_ps(
                _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(u_2, xi_173), xi_11),
                              xi_176),
                xi_281),
            _mm256_set_ps(xi_216, xi_216, xi_216, xi_216, xi_216, xi_216,
                          xi_216, xi_216));
        const __m256 xi_257 = _mm256_add_ps(
            _mm256_mul_ps(xi_255, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)),
            _mm256_mul_ps(xi_256, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)));
        const __m256 xi_262 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_mul_ps(xi_258,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_260,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0))),
                    xi_213),
                xi_259),
            xi_261);
        const __m256 xi_263 = _mm256_add_ps(xi_255, xi_256);
        const __m256 xi_264 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_mul_ps(xi_259,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_261,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0))),
                    xi_220),
                xi_258),
            xi_260);
        const __m256 xi_265 = _mm256_add_ps(
            _mm256_mul_ps(xi_239, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                -1.0, -1.0, -1.0)),
            xi_240);
        const __m256 xi_266 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_mul_ps(xi_234,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0)),
                        _mm256_mul_ps(xi_238,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0))),
                    xi_232),
                xi_236),
            xi_265);
        const __m256 xi_267 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(_mm256_add_ps(xi_199, xi_208), xi_250),
                        xi_251),
                    xi_252),
                xi_253),
            xi_265);
        const __m256 forceTerm_0 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_mul_ps(xi_26, _mm256_set_ps(-1.0, -1.0, -1.0,
                                                               -1.0, -1.0, -1.0,
                                                               -1.0, -1.0)),
                            _mm256_mul_ps(xi_26,
                                          _mm256_set_ps(xi_29, xi_29, xi_29,
                                                        xi_29, xi_29, xi_29,
                                                        xi_29, xi_29))),
                        _mm256_mul_ps(xi_27,
                                      _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                    -1.0, -1.0, -1.0, -1.0))),
                    _mm256_mul_ps(xi_27,
                                  _mm256_set_ps(xi_29, xi_29, xi_29, xi_29,
                                                xi_29, xi_29, xi_29, xi_29))),
                _mm256_mul_ps(xi_28, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0))),
            _mm256_mul_ps(xi_28, _mm256_set_ps(xi_29, xi_29, xi_29, xi_29,
                                               xi_29, xi_29, xi_29, xi_29)));
        const __m256 forceTerm_1 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_mul_ps(xi_32, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                xi_30),
            xi_42);
        const __m256 forceTerm_2 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_mul_ps(xi_30, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                xi_32),
            xi_42);
        const __m256 forceTerm_3 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_mul_ps(xi_43, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                xi_45),
            xi_49);
        const __m256 forceTerm_4 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_mul_ps(xi_45, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                xi_43),
            xi_49);
        const __m256 forceTerm_5 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_mul_ps(xi_52, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                xi_50),
            xi_53);
        const __m256 forceTerm_6 = _mm256_add_ps(
            _mm256_add_ps(
                _mm256_mul_ps(xi_50, _mm256_set_ps(-1.0, -1.0, -1.0, -1.0, -1.0,
                                                   -1.0, -1.0, -1.0)),
                xi_52),
            xi_53);
        const __m256 forceTerm_7 =
            _mm256_add_ps(_mm256_add_ps(xi_58, xi_67), xi_77);
        const __m256 forceTerm_8 =
            _mm256_add_ps(_mm256_add_ps(xi_67, xi_78), xi_79);
        const __m256 forceTerm_9 =
            _mm256_add_ps(_mm256_add_ps(xi_58, xi_79), xi_81);
        const __m256 forceTerm_10 =
            _mm256_add_ps(_mm256_add_ps(xi_77, xi_78), xi_81);
        const __m256 forceTerm_11 =
            _mm256_add_ps(_mm256_add_ps(xi_84, xi_86), xi_91);
        const __m256 forceTerm_12 =
            _mm256_add_ps(_mm256_add_ps(xi_84, xi_92), xi_93);
        const __m256 forceTerm_13 =
            _mm256_add_ps(_mm256_add_ps(xi_102, xi_58), xi_95);
        const __m256 forceTerm_14 =
            _mm256_add_ps(_mm256_add_ps(xi_103, xi_78), xi_95);
        const __m256 forceTerm_15 =
            _mm256_add_ps(_mm256_add_ps(xi_104, xi_86), xi_93);
        const __m256 forceTerm_16 =
            _mm256_add_ps(_mm256_add_ps(xi_104, xi_91), xi_92);
        const __m256 forceTerm_17 =
            _mm256_add_ps(_mm256_add_ps(xi_103, xi_105), xi_58);
        const __m256 forceTerm_18 =
            _mm256_add_ps(_mm256_add_ps(xi_102, xi_105), xi_78);
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
                                                    _mm256_set_ps(
                                                        0.142857142857143f,
                                                        0.142857142857143f,
                                                        0.142857142857143f,
                                                        0.142857142857143f,
                                                        0.142857142857143f,
                                                        0.142857142857143f,
                                                        0.142857142857143f,
                                                        0.142857142857143f)),
                                                _mm256_mul_ps(
                                                    xi_109,
                                                    _mm256_set_ps(0.2f, 0.2f,
                                                                  0.2f, 0.2f,
                                                                  0.2f, 0.2f,
                                                                  0.2f, 0.2f))),
                                            _mm256_mul_ps(
                                                xi_111,
                                                _mm256_set_ps(-1.0, -1.0, -1.0,
                                                              -1.0, -1.0, -1.0,
                                                              -1.0, -1.0))),
                                        _mm256_mul_ps(
                                            xi_112, _mm256_set_ps(
                                                        0.0857142857142857f,
                                                        0.0857142857142857f,
                                                        0.0857142857142857f,
                                                        0.0857142857142857f,
                                                        0.0857142857142857f,
                                                        0.0857142857142857f,
                                                        0.0857142857142857f,
                                                        0.0857142857142857f))),
                                    _mm256_mul_ps(
                                        xi_115,
                                        _mm256_set_ps(0.1f, 0.1f, 0.1f, 0.1f,
                                                      0.1f, 0.1f, 0.1f, 0.1f))),
                                _mm256_mul_ps(
                                    xi_118,
                                    _mm256_set_ps(0.0428571428571429f,
                                                  0.0428571428571429f,
                                                  0.0428571428571429f,
                                                  0.0428571428571429f,
                                                  0.0428571428571429f,
                                                  0.0428571428571429f,
                                                  0.0428571428571429f,
                                                  0.0428571428571429f))),
                            _mm256_mul_ps(xi_128,
                                          _mm256_set_ps(-0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f, -0.5f,
                                                        -0.5f, -0.5f))),
                        _mm256_mul_ps(
                            xi_132,
                            _mm256_set_ps(
                                0.0238095238095238f, 0.0238095238095238f,
                                0.0238095238095238f, 0.0238095238095238f,
                                0.0238095238095238f, 0.0238095238095238f,
                                0.0238095238095238f, 0.0238095238095238f))),
                    forceTerm_0),
                xi_274));
        _mm256_store_ps(
            &_data_pdfs_20_31_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_mul_ps(
                                        xi_136,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                                    _mm256_mul_ps(
                                        xi_146,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0))),
                                forceTerm_1),
                            xi_141),
                        xi_157),
                    xi_171),
                xi_288));
        _mm256_store_ps(
            &_data_pdfs_20_32_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_mul_ps(
                                        xi_141,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                                    forceTerm_2),
                                xi_136),
                            xi_146),
                        xi_171),
                    xi_172),
                xi_277));
        _mm256_store_ps(
            &_data_pdfs_20_33_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_mul_ps(
                                        xi_175,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                                    forceTerm_3),
                                xi_178),
                            xi_180),
                        xi_185),
                    xi_186),
                xi_273));
        _mm256_store_ps(
            &_data_pdfs_20_34_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_mul_ps(
                                        xi_178,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                                    _mm256_mul_ps(
                                        xi_180,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0))),
                                forceTerm_4),
                            xi_175),
                        xi_186),
                    xi_187),
                xi_282));
        _mm256_store_ps(
            &_data_pdfs_20_35_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_mul_ps(
                                        xi_190,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                                    _mm256_mul_ps(
                                        xi_194,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0))),
                                forceTerm_5),
                            xi_192),
                        xi_199),
                    xi_200),
                xi_272));
        _mm256_store_ps(
            &_data_pdfs_20_36_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_add_ps(
                            _mm256_add_ps(
                                _mm256_add_ps(
                                    _mm256_mul_ps(
                                        xi_192,
                                        _mm256_set_ps(-1.0, -1.0, -1.0, -1.0,
                                                      -1.0, -1.0, -1.0, -1.0)),
                                    forceTerm_6),
                                xi_190),
                            xi_194),
                        xi_200),
                    xi_201),
                xi_270));
        _mm256_store_ps(
            &_data_pdfs_20_37_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(forceTerm_7, xi_210), xi_214),
                    xi_218),
                xi_278));
        _mm256_store_ps(
            &_data_pdfs_20_38_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(forceTerm_8, xi_210), xi_219),
                    xi_221),
                xi_275));
        _mm256_store_ps(
            &_data_pdfs_20_39_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(forceTerm_9, xi_214), xi_219),
                    xi_223),
                xi_279));
        _mm256_store_ps(
            &_data_pdfs_20_310_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(forceTerm_10, xi_218), xi_221),
                    xi_223),
                xi_285));
        _mm256_store_ps(
            &_data_pdfs_20_311_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(forceTerm_11, xi_230), xi_242),
                    xi_247),
                xi_280));
        _mm256_store_ps(
            &_data_pdfs_20_312_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(forceTerm_12, xi_242), xi_248),
                    xi_249),
                xi_284));
        _mm256_store_ps(
            &_data_pdfs_20_313_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(forceTerm_13, xi_254), xi_257),
                    xi_262),
                xi_289));
        _mm256_store_ps(
            &_data_pdfs_20_314_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(forceTerm_14, xi_254), xi_263),
                    xi_264),
                xi_283));
        _mm256_store_ps(
            &_data_pdfs_20_315_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(forceTerm_15, xi_247), xi_249),
                    xi_266),
                xi_268));
        _mm256_store_ps(
            &_data_pdfs_20_316_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(forceTerm_16, xi_230), xi_248),
                    xi_266),
                xi_286));
        _mm256_store_ps(
            &_data_pdfs_20_317_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(forceTerm_17, xi_262), xi_263),
                    xi_267),
                xi_287));
        _mm256_store_ps(
            &_data_pdfs_20_318_10[ctr_0],
            _mm256_add_ps(
                _mm256_add_ps(
                    _mm256_add_ps(_mm256_add_ps(forceTerm_18, xi_257), xi_264),
                    xi_267),
                xi_281));
      }
    }
  }
}
} // namespace internal_48c9ee502281a70505dce0378c55abd5

void CollideSweepSinglePrecisionThermalizedAVX::run(IBlock *block) {
  auto force = block->getData<field::GhostLayerField<float, 3>>(forceID);
  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);

  auto &omega_even = this->omega_even_;
  auto block_offset_2 = this->block_offset_2_;
  auto &seed = this->seed_;
  auto &time_step = this->time_step_;
  auto block_offset_0 = this->block_offset_0_;
  auto &omega_shear = this->omega_shear_;
  auto block_offset_1 = this->block_offset_1_;
  auto &kT = this->kT_;
  auto &omega_odd = this->omega_odd_;
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
  internal_48c9ee502281a70505dce0378c55abd5::
      collidesweepsingleprecisionthermalizedavx_collidesweepsingleprecisionthermalizedavx(
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

  auto force = block->getData<field::GhostLayerField<float, 3>>(forceID);
  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);

  auto &omega_even = this->omega_even_;
  auto block_offset_2 = this->block_offset_2_;
  auto &seed = this->seed_;
  auto &time_step = this->time_step_;
  auto block_offset_0 = this->block_offset_0_;
  auto &omega_shear = this->omega_shear_;
  auto block_offset_1 = this->block_offset_1_;
  auto &kT = this->kT_;
  auto &omega_odd = this->omega_odd_;
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
  internal_48c9ee502281a70505dce0378c55abd5::
      collidesweepsingleprecisionthermalizedavx_collidesweepsingleprecisionthermalizedavx(
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
