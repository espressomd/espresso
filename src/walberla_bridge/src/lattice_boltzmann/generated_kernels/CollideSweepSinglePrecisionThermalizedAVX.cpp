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

// kernel generated with pystencils v1.2, lbmpy v1.2, lbmpy_walberla/pystencils_walberla from waLBerla commit 4d10e7f2358fc4a4f7e99195d0f67f0b759ecb6f

#include <cmath>

#include "CollideSweepSinglePrecisionThermalizedAVX.h"
#include "core/DataTypes.h"
#include "core/Macros.h"

#include "philox_rand.h"

#include <immintrin.h>

#define FUNC_PREFIX

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
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

namespace internal_48c9ee502281a70505dce0378c55abd5 {
static FUNC_PREFIX void collidesweepsingleprecisionthermalizedavx_collidesweepsingleprecisionthermalizedavx(float *RESTRICT const _data_force, float *RESTRICT _data_pdfs, int64_t const _size_force_0, int64_t const _size_force_1, int64_t const _size_force_2, int64_t const _stride_force_1, int64_t const _stride_force_2, int64_t const _stride_force_3, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, uint32_t block_offset_0, uint32_t block_offset_1, uint32_t block_offset_2, float kT, float omega_bulk, float omega_even, float omega_odd, float omega_shear, uint32_t seed, uint32_t time_step) {
  const float xi_28 = omega_bulk * 0.5f;
  const float xi_55 = omega_shear * 0.041666666666666664f;
  const float xi_60 = omega_bulk * 0.041666666666666664f;
  const float xi_71 = omega_shear * 0.125f;
  const float xi_109 = 2.4494897427831779f;
  const float xi_134 = omega_odd * 0.25f;
  const float xi_145 = omega_odd * 0.083333333333333329f;
  const float xi_198 = omega_shear * 0.25f;
  const float xi_211 = omega_odd * 0.041666666666666664f;
  const float xi_213 = omega_odd * 0.125f;
  const float rr_0 = 0.0f;
  const float xi_53 = rr_0 * 0.041666666666666664f;
  for (int64_t ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1) {
    float *RESTRICT _data_pdfs_20_36 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_32 = _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    float *RESTRICT _data_pdfs_20_31 = _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_32 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_311 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_318 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_313 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    float *RESTRICT _data_pdfs_20_34 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_33 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_317 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    float *RESTRICT _data_pdfs_20_35 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_314 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_312 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_316 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_38 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_315 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_31 = _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    float *RESTRICT _data_pdfs_20_310 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_39 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_37 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1) {
      float *RESTRICT _data_pdfs_20_36_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_36;
      float *RESTRICT _data_force_20_32_10 = _stride_force_1 * ctr_1 + _data_force_20_32;
      float *RESTRICT _data_pdfs_20_31_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_31;
      float *RESTRICT _data_pdfs_20_32_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_32;
      float *RESTRICT _data_pdfs_20_311_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_311;
      float *RESTRICT _data_pdfs_20_318_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_318;
      float *RESTRICT _data_pdfs_20_313_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_313;
      float *RESTRICT _data_pdfs_20_30_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      float *RESTRICT _data_pdfs_20_34_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      float *RESTRICT _data_pdfs_20_33_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      float *RESTRICT _data_pdfs_20_317_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_317;
      float *RESTRICT _data_force_20_30_10 = _stride_force_1 * ctr_1 + _data_force_20_30;
      float *RESTRICT _data_pdfs_20_35_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_35;
      float *RESTRICT _data_pdfs_20_314_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_314;
      float *RESTRICT _data_pdfs_20_312_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_312;
      float *RESTRICT _data_pdfs_20_316_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_316;
      float *RESTRICT _data_pdfs_20_38_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_38;
      float *RESTRICT _data_pdfs_20_315_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_315;
      float *RESTRICT _data_force_20_31_10 = _stride_force_1 * ctr_1 + _data_force_20_31;
      float *RESTRICT _data_pdfs_20_310_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_310;
      float *RESTRICT _data_pdfs_20_39_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_39;
      float *RESTRICT _data_pdfs_20_37_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_37;
      {
        for (int64_t ctr_0 = 0; ctr_0 < (int64_t)((_size_force_0) / (8)) * (8); ctr_0 += 8) {
          const __m256 xi_244 = _mm256_load_ps(&_data_pdfs_20_36_10[ctr_0]);
          const __m256 xi_245 = _mm256_load_ps(&_data_force_20_32_10[ctr_0]);
          const __m256 xi_246 = _mm256_load_ps(&_data_pdfs_20_31_10[ctr_0]);
          const __m256 xi_247 = _mm256_load_ps(&_data_pdfs_20_32_10[ctr_0]);
          const __m256 xi_248 = _mm256_load_ps(&_data_pdfs_20_311_10[ctr_0]);
          const __m256 xi_249 = _mm256_load_ps(&_data_pdfs_20_318_10[ctr_0]);
          const __m256 xi_250 = _mm256_load_ps(&_data_pdfs_20_313_10[ctr_0]);
          const __m256 xi_251 = _mm256_load_ps(&_data_pdfs_20_30_10[ctr_0]);
          const __m256 xi_252 = _mm256_load_ps(&_data_pdfs_20_34_10[ctr_0]);
          const __m256 xi_253 = _mm256_load_ps(&_data_pdfs_20_33_10[ctr_0]);
          const __m256 xi_254 = _mm256_load_ps(&_data_pdfs_20_317_10[ctr_0]);
          const __m256 xi_255 = _mm256_load_ps(&_data_force_20_30_10[ctr_0]);
          const __m256 xi_256 = _mm256_load_ps(&_data_pdfs_20_35_10[ctr_0]);
          const __m256 xi_257 = _mm256_load_ps(&_data_pdfs_20_314_10[ctr_0]);
          const __m256 xi_258 = _mm256_load_ps(&_data_pdfs_20_312_10[ctr_0]);
          const __m256 xi_259 = _mm256_load_ps(&_data_pdfs_20_316_10[ctr_0]);
          const __m256 xi_260 = _mm256_load_ps(&_data_pdfs_20_38_10[ctr_0]);
          const __m256 xi_261 = _mm256_load_ps(&_data_pdfs_20_315_10[ctr_0]);
          const __m256 xi_262 = _mm256_load_ps(&_data_force_20_31_10[ctr_0]);
          const __m256 xi_263 = _mm256_load_ps(&_data_pdfs_20_310_10[ctr_0]);
          const __m256 xi_264 = _mm256_load_ps(&_data_pdfs_20_39_10[ctr_0]);
          const __m256 xi_265 = _mm256_load_ps(&_data_pdfs_20_37_10[ctr_0]);

          __m256 random_3_0{};
          __m256 random_3_1{};
          __m256 random_3_2{};
          __m256 random_3_3{};
          if (kT > 0.) {
            philox_float4(time_step, _mm256_add_epi32(_mm256_add_epi32(_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0), _mm256_set_epi32(ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0)), _mm256_set_epi32(((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)))), block_offset_1 + ctr_1, block_offset_2 + ctr_2, 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);
          }

          __m256 random_2_0{};
          __m256 random_2_1{};
          __m256 random_2_2{};
          __m256 random_2_3{};
          if (kT > 0.) {
            philox_float4(time_step, _mm256_add_epi32(_mm256_add_epi32(_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0), _mm256_set_epi32(ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0)), _mm256_set_epi32(((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)))), block_offset_1 + ctr_1, block_offset_2 + ctr_2, 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);
          }

          __m256 random_1_0{};
          __m256 random_1_1{};
          __m256 random_1_2{};
          __m256 random_1_3{};
          if (kT > 0.) {
            philox_float4(time_step, _mm256_add_epi32(_mm256_add_epi32(_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0), _mm256_set_epi32(ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0)), _mm256_set_epi32(((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)))), block_offset_1 + ctr_1, block_offset_2 + ctr_2, 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);
          }

          __m256 random_0_0{};
          __m256 random_0_1{};
          __m256 random_0_2{};
          __m256 random_0_3{};
          if (kT > 0.) {
            philox_float4(time_step, _mm256_add_epi32(_mm256_add_epi32(_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0), _mm256_set_epi32(ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0)), _mm256_set_epi32(((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)))), block_offset_1 + ctr_1, block_offset_2 + ctr_2, 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);
          }
          const __m256 xi_2 = _mm256_add_ps(xi_249, xi_257);
          const __m256 xi_3 = _mm256_add_ps(xi_2, xi_252);
          const __m256 xi_4 = _mm256_add_ps(_mm256_add_ps(xi_246, xi_248), xi_261);
          const __m256 xi_5 = _mm256_add_ps(xi_256, xi_258);
          const __m256 xi_6 = _mm256_add_ps(xi_244, xi_254);
          const __m256 xi_8 = _mm256_mul_ps(xi_264, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f));
          const __m256 xi_9 = _mm256_mul_ps(xi_265, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f));
          const __m256 xi_10 = _mm256_mul_ps(xi_254, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f));
          const __m256 xi_11 = _mm256_mul_ps(xi_250, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f));
          const __m256 xi_12 = _mm256_mul_ps(xi_253, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f));
          const __m256 xi_13 = _mm256_add_ps(_mm256_add_ps(xi_10, xi_11), xi_12);
          const __m256 xi_14 = _mm256_mul_ps(xi_247, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f));
          const __m256 xi_15 = _mm256_mul_ps(xi_263, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f));
          const __m256 xi_16 = _mm256_add_ps(xi_14, xi_15);
          const __m256 xi_17 = _mm256_mul_ps(xi_259, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f));
          const __m256 xi_18 = _mm256_mul_ps(xi_258, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f));
          const __m256 xi_19 = _mm256_add_ps(xi_17, xi_18);
          const __m256 xi_20 = _mm256_mul_ps(xi_249, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f));
          const __m256 xi_21 = _mm256_add_ps(xi_10, xi_20);
          const __m256 xi_22 = _mm256_mul_ps(xi_261, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f));
          const __m256 xi_23 = _mm256_mul_ps(xi_244, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f));
          const __m256 xi_24 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_17, xi_22), xi_23), xi_248);
          const __m256 xi_29 = _mm256_mul_ps(xi_262, _mm256_set_ps(0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f));
          const __m256 xi_30 = _mm256_mul_ps(xi_262, _mm256_set_ps(0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f));
          const __m256 xi_42 = _mm256_mul_ps(xi_255, _mm256_set_ps(0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f));
          const __m256 xi_43 = _mm256_mul_ps(xi_255, _mm256_set_ps(0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f));
          const __m256 xi_49 = _mm256_mul_ps(xi_245, _mm256_set_ps(0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f));
          const __m256 xi_50 = _mm256_mul_ps(xi_245, _mm256_set_ps(0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f));
          const __m256 xi_67 = _mm256_mul_ps(xi_262, _mm256_set_ps(0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f));
          const __m256 xi_72 = _mm256_mul_ps(xi_262, _mm256_set_ps(xi_71, xi_71, xi_71, xi_71, xi_71, xi_71, xi_71, xi_71));
          const __m256 xi_114 = _mm256_mul_ps(xi_251, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f));
          const __m256 xi_118 = _mm256_mul_ps(xi_248, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f));
          const __m256 xi_119 = _mm256_add_ps(xi_118, xi_18);
          const __m256 xi_120 = _mm256_add_ps(_mm256_mul_ps(xi_260, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_8);
          const __m256 xi_122 = _mm256_mul_ps(xi_257, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f));
          const __m256 xi_123 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_11, xi_122), xi_15), xi_21);
          const __m256 xi_125 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_248, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f)), _mm256_mul_ps(xi_258, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f))), _mm256_mul_ps(xi_259, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f))), _mm256_mul_ps(xi_261, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f)));
          const __m256 xi_126 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_252, _mm256_set_ps(5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f)), _mm256_mul_ps(xi_253, _mm256_set_ps(5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f))), xi_125);
          const __m256 xi_128 = _mm256_mul_ps(xi_250, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f));
          const __m256 xi_129 = _mm256_mul_ps(xi_257, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f));
          const __m256 xi_130 = _mm256_add_ps(_mm256_mul_ps(xi_249, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f)), _mm256_mul_ps(xi_254, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f)));
          const __m256 xi_132 = _mm256_add_ps(xi_118, xi_258);
          const __m256 xi_133 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_132, xi_14), xi_22), xi_246), xi_259);
          const __m256 xi_135 = _mm256_mul_ps(xi_133, _mm256_set_ps(xi_134, xi_134, xi_134, xi_134, xi_134, xi_134, xi_134, xi_134));
          const __m256 xi_136 = _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_2_3);
          const __m256 xi_141 = _mm256_mul_ps(xi_265, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f));
          const __m256 xi_142 = _mm256_mul_ps(xi_263, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f));
          const __m256 xi_143 = _mm256_add_ps(_mm256_mul_ps(xi_264, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f)), _mm256_mul_ps(xi_260, _mm256_set_ps(-2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f)));
          const __m256 xi_144 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_141, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_14), xi_142), xi_143), xi_19), xi_4);
          const __m256 xi_146 = _mm256_mul_ps(xi_144, _mm256_set_ps(xi_145, xi_145, xi_145, xi_145, xi_145, xi_145, xi_145, xi_145));
          const __m256 xi_147 = _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_1_2);
          const __m256 xi_152 = _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_0_1);
          const __m256 xi_166 = _mm256_add_ps(xi_122, xi_250);
          const __m256 xi_167 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_12, xi_166), xi_20), xi_252), xi_254);
          const __m256 xi_168 = _mm256_mul_ps(xi_167, _mm256_set_ps(xi_134, xi_134, xi_134, xi_134, xi_134, xi_134, xi_134, xi_134));
          const __m256 xi_169 = _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_2_1);
          const __m256 xi_171 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_142, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_13), xi_141), xi_143), xi_3);
          const __m256 xi_172 = _mm256_mul_ps(xi_171, _mm256_set_ps(xi_145, xi_145, xi_145, xi_145, xi_145, xi_145, xi_145, xi_145));
          const __m256 xi_173 = _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_2_0);
          const __m256 xi_178 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_119, xi_23), xi_256), xi_259), xi_261);
          const __m256 xi_179 = _mm256_mul_ps(xi_178, _mm256_set_ps(xi_134, xi_134, xi_134, xi_134, xi_134, xi_134, xi_134, xi_134));
          const __m256 xi_180 = _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_2_2);
          const __m256 xi_182 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_128, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_129, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_130), xi_24), xi_5);
          const __m256 xi_183 = _mm256_mul_ps(xi_182, _mm256_set_ps(xi_145, xi_145, xi_145, xi_145, xi_145, xi_145, xi_145, xi_145));
          const __m256 xi_184 = _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_1_3);
          const __m256 xi_212 = _mm256_mul_ps(xi_182, _mm256_set_ps(xi_211, xi_211, xi_211, xi_211, xi_211, xi_211, xi_211, xi_211));
          const __m256 xi_214 = _mm256_mul_ps(xi_178, _mm256_set_ps(xi_213, xi_213, xi_213, xi_213, xi_213, xi_213, xi_213, xi_213));
          const __m256 xi_220 = _mm256_mul_ps(xi_144, _mm256_set_ps(xi_211, xi_211, xi_211, xi_211, xi_211, xi_211, xi_211, xi_211));
          const __m256 xi_221 = _mm256_mul_ps(xi_133, _mm256_set_ps(xi_213, xi_213, xi_213, xi_213, xi_213, xi_213, xi_213, xi_213));
          const __m256 xi_235 = _mm256_mul_ps(xi_167, _mm256_set_ps(xi_213, xi_213, xi_213, xi_213, xi_213, xi_213, xi_213, xi_213));
          const __m256 xi_236 = _mm256_mul_ps(xi_171, _mm256_set_ps(xi_211, xi_211, xi_211, xi_211, xi_211, xi_211, xi_211, xi_211));
          const __m256 xi_31 = _mm256_mul_ps(xi_30, _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0, rr_0, rr_0, rr_0));
          const __m256 xi_44 = _mm256_mul_ps(xi_43, _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0, rr_0, rr_0, rr_0));
          const __m256 xi_51 = _mm256_mul_ps(xi_50, _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0, rr_0, rr_0, rr_0));
          const __m256 xi_54 = _mm256_mul_ps(xi_255, _mm256_set_ps(xi_53, xi_53, xi_53, xi_53, xi_53, xi_53, xi_53, xi_53));
          const __m256 xi_59 = _mm256_mul_ps(xi_262, _mm256_set_ps(xi_53, xi_53, xi_53, xi_53, xi_53, xi_53, xi_53, xi_53));
          const __m256 xi_81 = _mm256_mul_ps(xi_245, _mm256_set_ps(xi_53, xi_53, xi_53, xi_53, xi_53, xi_53, xi_53, xi_53));
          const __m256 vel0Term = _mm256_add_ps(_mm256_add_ps(xi_260, xi_263), xi_3);
          const __m256 vel1Term = _mm256_add_ps(xi_265, xi_4);
          const __m256 vel2Term = _mm256_add_ps(xi_250, xi_5);
          const __m256 rho = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(vel0Term, vel1Term), vel2Term), xi_247), xi_251), xi_253), xi_259), xi_264), xi_6);
          const __m256 xi_105 = _mm256_mul_ps(rho, _mm256_set_ps(kT, kT, kT, kT, kT, kT, kT, kT));
          const __m256 xi_106 = _mm256_sqrt_ps(_mm256_mul_ps(xi_105, _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_mul_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_set_ps(omega_even, omega_even, omega_even, omega_even, omega_even, omega_even, omega_even, omega_even)), _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f)), _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_set_ps(omega_even, omega_even, omega_even, omega_even, omega_even, omega_even, omega_even, omega_even)), _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f)))), _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f))));
          const __m256 xi_107 = _mm256_mul_ps(_mm256_mul_ps(xi_106, _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_3_0)), _mm256_set_ps(3.7416573867739413f, 3.7416573867739413f, 3.7416573867739413f, 3.7416573867739413f, 3.7416573867739413f, 3.7416573867739413f, 3.7416573867739413f, 3.7416573867739413f));
          const __m256 xi_108 = _mm256_mul_ps(_mm256_mul_ps(xi_106, _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_3_2)), _mm256_set_ps(5.4772255750516612f, 5.4772255750516612f, 5.4772255750516612f, 5.4772255750516612f, 5.4772255750516612f, 5.4772255750516612f, 5.4772255750516612f, 5.4772255750516612f));
          const __m256 xi_110 = _mm256_mul_ps(_mm256_mul_ps(_mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_1_1), _mm256_set_ps(xi_109, xi_109, xi_109, xi_109, xi_109, xi_109, xi_109, xi_109)), _mm256_sqrt_ps(_mm256_mul_ps(xi_105, _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_mul_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk)), _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f)), _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk)), _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f)))), _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f)))));
          const __m256 xi_111 = _mm256_mul_ps(_mm256_mul_ps(xi_106, _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_3_1)), _mm256_set_ps(8.3666002653407556f, 8.3666002653407556f, 8.3666002653407556f, 8.3666002653407556f, 8.3666002653407556f, 8.3666002653407556f, 8.3666002653407556f, 8.3666002653407556f));
          const __m256 xi_137 = _mm256_sqrt_ps(_mm256_mul_ps(xi_105, _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_mul_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_set_ps(omega_odd, omega_odd, omega_odd, omega_odd, omega_odd, omega_odd, omega_odd, omega_odd)), _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f)), _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_set_ps(omega_odd, omega_odd, omega_odd, omega_odd, omega_odd, omega_odd, omega_odd, omega_odd)), _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f)))), _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f))));
          const __m256 xi_138 = _mm256_mul_ps(xi_137, _mm256_set_ps(1.4142135623730951f, 1.4142135623730951f, 1.4142135623730951f, 1.4142135623730951f, 1.4142135623730951f, 1.4142135623730951f, 1.4142135623730951f, 1.4142135623730951f));
          const __m256 xi_139 = _mm256_mul_ps(xi_138, _mm256_set_ps(0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f));
          const __m256 xi_140 = _mm256_mul_ps(xi_136, xi_139);
          const __m256 xi_148 = _mm256_mul_ps(xi_137, _mm256_set_ps(xi_109, xi_109, xi_109, xi_109, xi_109, xi_109, xi_109, xi_109));
          const __m256 xi_149 = _mm256_mul_ps(xi_148, _mm256_set_ps(0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f));
          const __m256 xi_150 = _mm256_mul_ps(xi_147, xi_149);
          const __m256 xi_151 = _mm256_add_ps(_mm256_mul_ps(xi_146, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_150, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)));
          const __m256 xi_153 = _mm256_sqrt_ps(_mm256_mul_ps(xi_105, _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_mul_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear)), _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f)), _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear)), _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f)))), _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f))));
          const __m256 xi_154 = _mm256_mul_ps(xi_153, _mm256_set_ps(0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f));
          const __m256 xi_155 = _mm256_mul_ps(xi_152, xi_154);
          const __m256 xi_161 = _mm256_mul_ps(_mm256_mul_ps(xi_153, _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_0_0)), _mm256_set_ps(1.7320508075688772f, 1.7320508075688772f, 1.7320508075688772f, 1.7320508075688772f, 1.7320508075688772f, 1.7320508075688772f, 1.7320508075688772f, 1.7320508075688772f));
          const __m256 xi_165 = _mm256_add_ps(xi_146, xi_150);
          const __m256 xi_170 = _mm256_mul_ps(xi_139, xi_169);
          const __m256 xi_174 = _mm256_mul_ps(xi_149, xi_173);
          const __m256 xi_175 = _mm256_add_ps(xi_172, xi_174);
          const __m256 xi_177 = _mm256_add_ps(_mm256_mul_ps(xi_172, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_174, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)));
          const __m256 xi_181 = _mm256_mul_ps(xi_139, xi_180);
          const __m256 xi_185 = _mm256_mul_ps(xi_149, xi_184);
          const __m256 xi_186 = _mm256_add_ps(_mm256_mul_ps(xi_183, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_185, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)));
          const __m256 xi_188 = _mm256_add_ps(xi_183, xi_185);
          const __m256 xi_189 = _mm256_mul_ps(_mm256_mul_ps(xi_152, xi_153), _mm256_set_ps(0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f));
          const __m256 xi_192 = _mm256_mul_ps(xi_107, _mm256_set_ps(0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f));
          const __m256 xi_196 = _mm256_mul_ps(xi_154, _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_0_2));
          const __m256 xi_203 = _mm256_mul_ps(xi_154, _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_1_0));
          const __m256 xi_207 = _mm256_mul_ps(xi_111, _mm256_set_ps(-0.014285714285714285f, -0.014285714285714285f, -0.014285714285714285f, -0.014285714285714285f, -0.014285714285714285f, -0.014285714285714285f, -0.014285714285714285f, -0.014285714285714285f));
          const __m256 xi_208 = _mm256_mul_ps(xi_108, _mm256_set_ps(0.050000000000000003f, 0.050000000000000003f, 0.050000000000000003f, 0.050000000000000003f, 0.050000000000000003f, 0.050000000000000003f, 0.050000000000000003f, 0.050000000000000003f));
          const __m256 xi_215 = _mm256_mul_ps(xi_148, _mm256_set_ps(0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f));
          const __m256 xi_216 = _mm256_mul_ps(xi_184, xi_215);
          const __m256 xi_217 = _mm256_mul_ps(xi_138, _mm256_set_ps(0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f));
          const __m256 xi_218 = _mm256_mul_ps(xi_180, xi_217);
          const __m256 xi_219 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_212, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_216, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_214), xi_218);
          const __m256 xi_222 = _mm256_mul_ps(xi_147, xi_215);
          const __m256 xi_223 = _mm256_mul_ps(xi_136, xi_217);
          const __m256 xi_224 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_220, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_222, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_221), xi_223);
          const __m256 xi_225 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_221, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_223, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_220), xi_222);
          const __m256 xi_227 = _mm256_mul_ps(xi_189, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f));
          const __m256 xi_230 = _mm256_mul_ps(xi_111, _mm256_set_ps(0.035714285714285712f, 0.035714285714285712f, 0.035714285714285712f, 0.035714285714285712f, 0.035714285714285712f, 0.035714285714285712f, 0.035714285714285712f, 0.035714285714285712f));
          const __m256 xi_232 = _mm256_mul_ps(xi_154, _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_0_3));
          const __m256 xi_237 = _mm256_mul_ps(xi_169, xi_217);
          const __m256 xi_238 = _mm256_mul_ps(xi_173, xi_215);
          const __m256 xi_239 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_235, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_237, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_236), xi_238);
          const __m256 xi_241 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_236, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_238, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_235), xi_237);
          const __m256 xi_242 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_214, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_218, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_212), xi_216);
          const __m256 xi_0 = _mm256_div_ps(_mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f), rho);
          const __m256 xi_7 = _mm256_mul_ps(xi_0, _mm256_set_ps(0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f));
          const __m256 u_0 = _mm256_add_ps(_mm256_mul_ps(xi_0, _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(vel0Term, xi_13), xi_8), xi_9)), _mm256_mul_ps(xi_255, xi_7));
          const __m256 xi_25 = _mm256_mul_ps(u_0, xi_255);
          const __m256 xi_37 = _mm256_mul_ps(xi_25, _mm256_set_ps(0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f));
          const __m256 xi_38 = _mm256_mul_ps(xi_25, _mm256_set_ps(0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f));
          const __m256 xi_39 = _mm256_mul_ps(xi_38, _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear));
          const __m256 xi_40 = _mm256_add_ps(_mm256_mul_ps(xi_37, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_39);
          const __m256 xi_56 = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(xi_25, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_set_ps(xi_55, xi_55, xi_55, xi_55, xi_55, xi_55, xi_55, xi_55)), xi_37);
          const __m256 xi_57 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_43, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_54), xi_56);
          const __m256 xi_61 = _mm256_mul_ps(_mm256_mul_ps(xi_25, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_set_ps(xi_60, xi_60, xi_60, xi_60, xi_60, xi_60, xi_60, xi_60));
          const __m256 xi_68 = _mm256_mul_ps(u_0, xi_67);
          const __m256 xi_73 = _mm256_mul_ps(u_0, xi_72);
          const __m256 xi_77 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_54, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_43), xi_56);
          const __m256 xi_84 = _mm256_mul_ps(xi_38, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f));
          const __m256 xi_95 = _mm256_mul_ps(u_0, xi_245);
          const __m256 xi_96 = _mm256_mul_ps(xi_95, _mm256_set_ps(0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f));
          const __m256 xi_99 = _mm256_mul_ps(xi_95, _mm256_set_ps(xi_71, xi_71, xi_71, xi_71, xi_71, xi_71, xi_71, xi_71));
          const __m256 xi_113 = _mm256_mul_ps(rho, _mm256_mul_ps(u_0, u_0));
          const __m256 u_1 = _mm256_add_ps(_mm256_mul_ps(xi_0, _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(vel1Term, xi_16), xi_19), xi_260), xi_8)), _mm256_mul_ps(xi_262, xi_7));
          const __m256 xi_26 = _mm256_mul_ps(u_1, xi_262);
          const __m256 xi_32 = _mm256_mul_ps(xi_26, _mm256_set_ps(0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f));
          const __m256 xi_45 = _mm256_mul_ps(xi_26, _mm256_set_ps(0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f));
          const __m256 xi_46 = _mm256_mul_ps(xi_45, _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear));
          const __m256 xi_47 = _mm256_add_ps(_mm256_mul_ps(xi_32, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_46);
          const __m256 xi_62 = _mm256_mul_ps(_mm256_mul_ps(xi_26, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_set_ps(xi_60, xi_60, xi_60, xi_60, xi_60, xi_60, xi_60, xi_60));
          const __m256 xi_69 = _mm256_mul_ps(u_1, _mm256_set_ps(0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f));
          const __m256 xi_70 = _mm256_mul_ps(xi_255, xi_69);
          const __m256 xi_74 = _mm256_mul_ps(u_1, _mm256_set_ps(xi_71, xi_71, xi_71, xi_71, xi_71, xi_71, xi_71, xi_71));
          const __m256 xi_75 = _mm256_mul_ps(xi_255, xi_74);
          const __m256 xi_76 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_68, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_70, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_73), xi_75);
          const __m256 xi_78 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_73, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_75, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_68), xi_70);
          const __m256 xi_86 = _mm256_mul_ps(xi_245, xi_69);
          const __m256 xi_88 = _mm256_mul_ps(xi_245, xi_74);
          const __m256 xi_93 = _mm256_mul_ps(xi_45, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f));
          const __m256 xi_112 = _mm256_mul_ps(rho, _mm256_mul_ps(u_1, u_1));
          const __m256 xi_121 = _mm256_add_ps(_mm256_add_ps(xi_112, xi_120), xi_9);
          const __m256 xi_197 = _mm256_mul_ps(rho, u_1);
          const __m256 xi_199 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(u_0, xi_197), xi_120), xi_263), xi_265), _mm256_set_ps(xi_198, xi_198, xi_198, xi_198, xi_198, xi_198, xi_198, xi_198));
          const __m256 xi_200 = _mm256_add_ps(_mm256_mul_ps(xi_196, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_199, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)));
          const __m256 xi_201 = _mm256_add_ps(xi_196, xi_199);
          const __m256 u_2 = _mm256_add_ps(_mm256_mul_ps(xi_0, _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(vel2Term, xi_21), xi_24), xi_257)), _mm256_mul_ps(xi_245, xi_7));
          const __m256 xi_27 = _mm256_mul_ps(u_2, xi_245);
          const __m256 xi_33 = _mm256_mul_ps(xi_27, _mm256_set_ps(0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f));
          const __m256 xi_34 = _mm256_mul_ps(xi_27, _mm256_set_ps(0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f));
          const __m256 xi_35 = _mm256_mul_ps(xi_34, _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear));
          const __m256 xi_36 = _mm256_add_ps(_mm256_mul_ps(xi_33, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_35);
          const __m256 xi_41 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_26, _mm256_set_ps(0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f)), _mm256_mul_ps(_mm256_mul_ps(xi_32, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear))), xi_36), xi_40);
          const __m256 xi_48 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_25, _mm256_set_ps(0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f)), _mm256_mul_ps(_mm256_mul_ps(xi_37, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear))), xi_36), xi_47);
          const __m256 xi_52 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_27, _mm256_set_ps(0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f)), _mm256_mul_ps(_mm256_mul_ps(xi_33, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear))), xi_40), xi_47);
          const __m256 xi_58 = _mm256_mul_ps(xi_34, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f));
          const __m256 xi_63 = _mm256_mul_ps(_mm256_mul_ps(xi_27, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_set_ps(xi_60, xi_60, xi_60, xi_60, xi_60, xi_60, xi_60, xi_60));
          const __m256 xi_64 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(xi_26, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_set_ps(xi_55, xi_55, xi_55, xi_55, xi_55, xi_55, xi_55, xi_55)), xi_32), xi_61), xi_62), xi_63);
          const __m256 xi_65 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_59, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_30), xi_64);
          const __m256 xi_66 = _mm256_add_ps(_mm256_add_ps(xi_35, xi_58), xi_65);
          const __m256 xi_79 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_30, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_59), xi_64);
          const __m256 xi_80 = _mm256_add_ps(_mm256_add_ps(xi_35, xi_58), xi_79);
          const __m256 xi_82 = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(xi_27, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_set_ps(xi_55, xi_55, xi_55, xi_55, xi_55, xi_55, xi_55, xi_55)), xi_33);
          const __m256 xi_83 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_81, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_50), xi_82);
          const __m256 xi_85 = _mm256_add_ps(_mm256_add_ps(xi_39, xi_65), xi_84);
          const __m256 xi_87 = _mm256_mul_ps(u_2, xi_67);
          const __m256 xi_89 = _mm256_mul_ps(u_2, xi_72);
          const __m256 xi_90 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_88, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_89, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_86), xi_87);
          const __m256 xi_91 = _mm256_add_ps(_mm256_add_ps(xi_39, xi_79), xi_84);
          const __m256 xi_92 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_86, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_87, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_88), xi_89);
          const __m256 xi_94 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_46, xi_61), xi_62), xi_63), xi_83), xi_93);
          const __m256 xi_97 = _mm256_mul_ps(u_2, xi_255);
          const __m256 xi_98 = _mm256_mul_ps(xi_97, _mm256_set_ps(0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f));
          const __m256 xi_100 = _mm256_mul_ps(xi_97, _mm256_set_ps(xi_71, xi_71, xi_71, xi_71, xi_71, xi_71, xi_71, xi_71));
          const __m256 xi_101 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_96, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_98, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_100), xi_99);
          const __m256 xi_102 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_100, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_99, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_96), xi_98);
          const __m256 xi_103 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_50, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_81), xi_82);
          const __m256 xi_104 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_103, xi_46), xi_61), xi_62), xi_63), xi_93);
          const __m256 xi_115 = _mm256_mul_ps(rho, _mm256_mul_ps(u_2, u_2));
          const __m256 xi_116 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_244, _mm256_set_ps(3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f)), _mm256_mul_ps(xi_256, _mm256_set_ps(3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f))), _mm256_mul_ps(xi_115, _mm256_set_ps(0.66666666666666663f, 0.66666666666666663f, 0.66666666666666663f, 0.66666666666666663f, 0.66666666666666663f, 0.66666666666666663f, 0.66666666666666663f, 0.66666666666666663f))), xi_114);
          const __m256 xi_117 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_246, _mm256_set_ps(3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f)), _mm256_mul_ps(xi_247, _mm256_set_ps(3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f))), _mm256_mul_ps(xi_112, _mm256_set_ps(0.66666666666666663f, 0.66666666666666663f, 0.66666666666666663f, 0.66666666666666663f, 0.66666666666666663f, 0.66666666666666663f, 0.66666666666666663f, 0.66666666666666663f))), _mm256_mul_ps(xi_113, _mm256_set_ps(1.6666666666666667f, 1.6666666666666667f, 1.6666666666666667f, 1.6666666666666667f, 1.6666666666666667f, 1.6666666666666667f, 1.6666666666666667f, 1.6666666666666667f))), _mm256_mul_ps(xi_248, _mm256_set_ps(-3.0f, -3.0f, -3.0f, -3.0f, -3.0f, -3.0f, -3.0f, -3.0f))), _mm256_mul_ps(xi_258, _mm256_set_ps(-3.0f, -3.0f, -3.0f, -3.0f, -3.0f, -3.0f, -3.0f, -3.0f))), _mm256_mul_ps(xi_259, _mm256_set_ps(-3.0f, -3.0f, -3.0f, -3.0f, -3.0f, -3.0f, -3.0f, -3.0f))), _mm256_mul_ps(xi_261, _mm256_set_ps(-3.0f, -3.0f, -3.0f, -3.0f, -3.0f, -3.0f, -3.0f, -3.0f))), xi_116), _mm256_set_ps(omega_even, omega_even, omega_even, omega_even, omega_even, omega_even, omega_even, omega_even));
          const __m256 xi_124 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_113, xi_115), xi_119), xi_121), xi_123), xi_17), xi_22), xi_251), _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk));
          const __m256 xi_127 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_112, _mm256_set_ps(2.3333333333333335f, 2.3333333333333335f, 2.3333333333333335f, 2.3333333333333335f, 2.3333333333333335f, 2.3333333333333335f, 2.3333333333333335f, 2.3333333333333335f)), _mm256_mul_ps(xi_246, _mm256_set_ps(-2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f))), _mm256_mul_ps(xi_247, _mm256_set_ps(-2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f))), _mm256_mul_ps(xi_249, _mm256_set_ps(-5.0f, -5.0f, -5.0f, -5.0f, -5.0f, -5.0f, -5.0f, -5.0f))), _mm256_mul_ps(xi_250, _mm256_set_ps(-5.0f, -5.0f, -5.0f, -5.0f, -5.0f, -5.0f, -5.0f, -5.0f))), _mm256_mul_ps(xi_254, _mm256_set_ps(-5.0f, -5.0f, -5.0f, -5.0f, -5.0f, -5.0f, -5.0f, -5.0f))), _mm256_mul_ps(xi_257, _mm256_set_ps(-5.0f, -5.0f, -5.0f, -5.0f, -5.0f, -5.0f, -5.0f, -5.0f))), xi_116), xi_126), _mm256_set_ps(omega_even, omega_even, omega_even, omega_even, omega_even, omega_even, omega_even, omega_even));
          const __m256 xi_131 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_115, _mm256_set_ps(3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f)), _mm256_mul_ps(xi_246, _mm256_set_ps(5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f))), _mm256_mul_ps(xi_247, _mm256_set_ps(5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f))), _mm256_mul_ps(xi_244, _mm256_set_ps(-4.0f, -4.0f, -4.0f, -4.0f, -4.0f, -4.0f, -4.0f, -4.0f))), _mm256_mul_ps(xi_256, _mm256_set_ps(-4.0f, -4.0f, -4.0f, -4.0f, -4.0f, -4.0f, -4.0f, -4.0f))), _mm256_mul_ps(xi_260, _mm256_set_ps(-7.0f, -7.0f, -7.0f, -7.0f, -7.0f, -7.0f, -7.0f, -7.0f))), _mm256_mul_ps(xi_263, _mm256_set_ps(-7.0f, -7.0f, -7.0f, -7.0f, -7.0f, -7.0f, -7.0f, -7.0f))), _mm256_mul_ps(xi_264, _mm256_set_ps(-7.0f, -7.0f, -7.0f, -7.0f, -7.0f, -7.0f, -7.0f, -7.0f))), _mm256_mul_ps(xi_265, _mm256_set_ps(-7.0f, -7.0f, -7.0f, -7.0f, -7.0f, -7.0f, -7.0f, -7.0f))), xi_114), xi_126), xi_128), xi_129), xi_130), _mm256_set_ps(omega_even, omega_even, omega_even, omega_even, omega_even, omega_even, omega_even, omega_even));
          const __m256 xi_156 = _mm256_add_ps(_mm256_mul_ps(xi_115, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_256);
          const __m256 xi_157 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_246, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_121), xi_156), xi_16), xi_2), xi_250), xi_6), _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear));
          const __m256 xi_158 = _mm256_mul_ps(xi_157, _mm256_set_ps(0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f));
          const __m256 xi_159 = _mm256_add_ps(_mm256_mul_ps(xi_131, _mm256_set_ps(-0.01984126984126984f, -0.01984126984126984f, -0.01984126984126984f, -0.01984126984126984f, -0.01984126984126984f, -0.01984126984126984f, -0.01984126984126984f, -0.01984126984126984f)), _mm256_mul_ps(xi_107, _mm256_set_ps(-0.11904761904761904f, -0.11904761904761904f, -0.11904761904761904f, -0.11904761904761904f, -0.11904761904761904f, -0.11904761904761904f, -0.11904761904761904f, -0.11904761904761904f)));
          const __m256 xi_160 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_112, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_113, _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f))), _mm256_mul_ps(xi_252, _mm256_set_ps(-2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f))), _mm256_mul_ps(xi_253, _mm256_set_ps(-2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f))), xi_120), xi_123), xi_125), xi_156), xi_244), xi_246), xi_247), xi_9), _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear));
          const __m256 xi_162 = _mm256_add_ps(_mm256_mul_ps(xi_160, _mm256_set_ps(-0.041666666666666664f, -0.041666666666666664f, -0.041666666666666664f, -0.041666666666666664f, -0.041666666666666664f, -0.041666666666666664f, -0.041666666666666664f, -0.041666666666666664f)), _mm256_mul_ps(xi_161, _mm256_set_ps(-0.16666666666666666f, -0.16666666666666666f, -0.16666666666666666f, -0.16666666666666666f, -0.16666666666666666f, -0.16666666666666666f, -0.16666666666666666f, -0.16666666666666666f)));
          const __m256 xi_163 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_117, _mm256_set_ps(-0.050000000000000003f, -0.050000000000000003f, -0.050000000000000003f, -0.050000000000000003f, -0.050000000000000003f, -0.050000000000000003f, -0.050000000000000003f, -0.050000000000000003f)), _mm256_mul_ps(xi_108, _mm256_set_ps(-0.10000000000000001f, -0.10000000000000001f, -0.10000000000000001f, -0.10000000000000001f, -0.10000000000000001f, -0.10000000000000001f, -0.10000000000000001f, -0.10000000000000001f))), xi_162);
          const __m256 xi_164 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_127, _mm256_set_ps(0.014285714285714285f, 0.014285714285714285f, 0.014285714285714285f, 0.014285714285714285f, 0.014285714285714285f, 0.014285714285714285f, 0.014285714285714285f, 0.014285714285714285f)), _mm256_mul_ps(xi_111, _mm256_set_ps(0.028571428571428571f, 0.028571428571428571f, 0.028571428571428571f, 0.028571428571428571f, 0.028571428571428571f, 0.028571428571428571f, 0.028571428571428571f, 0.028571428571428571f))), xi_155), xi_158), xi_159), xi_163);
          const __m256 xi_176 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_160, _mm256_set_ps(0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f)), _mm256_mul_ps(xi_161, _mm256_set_ps(0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f))), _mm256_mul_ps(xi_127, _mm256_set_ps(-0.035714285714285712f, -0.035714285714285712f, -0.035714285714285712f, -0.035714285714285712f, -0.035714285714285712f, -0.035714285714285712f, -0.035714285714285712f, -0.035714285714285712f))), _mm256_mul_ps(xi_111, _mm256_set_ps(-0.071428571428571425f, -0.071428571428571425f, -0.071428571428571425f, -0.071428571428571425f, -0.071428571428571425f, -0.071428571428571425f, -0.071428571428571425f, -0.071428571428571425f))), xi_159);
          const __m256 xi_187 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_155, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_158, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_131, _mm256_set_ps(0.015873015873015872f, 0.015873015873015872f, 0.015873015873015872f, 0.015873015873015872f, 0.015873015873015872f, 0.015873015873015872f, 0.015873015873015872f, 0.015873015873015872f))), _mm256_mul_ps(xi_107, _mm256_set_ps(0.095238095238095233f, 0.095238095238095233f, 0.095238095238095233f, 0.095238095238095233f, 0.095238095238095233f, 0.095238095238095233f, 0.095238095238095233f, 0.095238095238095233f))), _mm256_mul_ps(xi_127, _mm256_set_ps(-0.021428571428571429f, -0.021428571428571429f, -0.021428571428571429f, -0.021428571428571429f, -0.021428571428571429f, -0.021428571428571429f, -0.021428571428571429f, -0.021428571428571429f))), _mm256_mul_ps(xi_111, _mm256_set_ps(-0.042857142857142858f, -0.042857142857142858f, -0.042857142857142858f, -0.042857142857142858f, -0.042857142857142858f, -0.042857142857142858f, -0.042857142857142858f, -0.042857142857142858f))), xi_163);
          const __m256 xi_190 = _mm256_mul_ps(xi_157, _mm256_set_ps(0.0625f, 0.0625f, 0.0625f, 0.0625f, 0.0625f, 0.0625f, 0.0625f, 0.0625f));
          const __m256 xi_191 = _mm256_mul_ps(xi_131, _mm256_set_ps(0.013888888888888888f, 0.013888888888888888f, 0.013888888888888888f, 0.013888888888888888f, 0.013888888888888888f, 0.013888888888888888f, 0.013888888888888888f, 0.013888888888888888f));
          const __m256 xi_193 = _mm256_add_ps(_mm256_mul_ps(xi_124, _mm256_set_ps(0.041666666666666664f, 0.041666666666666664f, 0.041666666666666664f, 0.041666666666666664f, 0.041666666666666664f, 0.041666666666666664f, 0.041666666666666664f, 0.041666666666666664f)), _mm256_mul_ps(xi_110, _mm256_set_ps(0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f)));
          const __m256 xi_194 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_160, _mm256_set_ps(0.020833333333333332f, 0.020833333333333332f, 0.020833333333333332f, 0.020833333333333332f, 0.020833333333333332f, 0.020833333333333332f, 0.020833333333333332f, 0.020833333333333332f)), _mm256_mul_ps(xi_161, _mm256_set_ps(0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f))), xi_193);
          const __m256 xi_195 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_165, xi_189), xi_190), xi_191), xi_192), xi_194);
          const __m256 xi_202 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_151, xi_189), xi_190), xi_191), xi_192), xi_194);
          const __m256 xi_204 = _mm256_mul_ps(xi_127, _mm256_set_ps(-0.0071428571428571426f, -0.0071428571428571426f, -0.0071428571428571426f, -0.0071428571428571426f, -0.0071428571428571426f, -0.0071428571428571426f, -0.0071428571428571426f, -0.0071428571428571426f));
          const __m256 xi_205 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(u_2, xi_197), xi_132), xi_17), xi_261), _mm256_set_ps(xi_198, xi_198, xi_198, xi_198, xi_198, xi_198, xi_198, xi_198));
          const __m256 xi_206 = _mm256_mul_ps(xi_117, _mm256_set_ps(0.025000000000000001f, 0.025000000000000001f, 0.025000000000000001f, 0.025000000000000001f, 0.025000000000000001f, 0.025000000000000001f, 0.025000000000000001f, 0.025000000000000001f));
          const __m256 xi_209 = _mm256_add_ps(_mm256_mul_ps(xi_131, _mm256_set_ps(-0.003968253968253968f, -0.003968253968253968f, -0.003968253968253968f, -0.003968253968253968f, -0.003968253968253968f, -0.003968253968253968f, -0.003968253968253968f, -0.003968253968253968f)), _mm256_mul_ps(xi_107, _mm256_set_ps(-0.023809523809523808f, -0.023809523809523808f, -0.023809523809523808f, -0.023809523809523808f, -0.023809523809523808f, -0.023809523809523808f, -0.023809523809523808f, -0.023809523809523808f)));
          const __m256 xi_210 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_162, xi_193), xi_203), xi_204), xi_205), xi_206), xi_207), xi_208), xi_209);
          const __m256 xi_226 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_203, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_205, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_162), xi_193), xi_204), xi_206), xi_207), xi_208), xi_209);
          const __m256 xi_228 = _mm256_mul_ps(xi_190, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f));
          const __m256 xi_229 = _mm256_mul_ps(xi_127, _mm256_set_ps(0.017857142857142856f, 0.017857142857142856f, 0.017857142857142856f, 0.017857142857142856f, 0.017857142857142856f, 0.017857142857142856f, 0.017857142857142856f, 0.017857142857142856f));
          const __m256 xi_231 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_188, xi_194), xi_209), xi_227), xi_228), xi_229), xi_230);
          const __m256 xi_233 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(rho, u_0), u_2), xi_10), xi_166), xi_249), _mm256_set_ps(xi_198, xi_198, xi_198, xi_198, xi_198, xi_198, xi_198, xi_198));
          const __m256 xi_234 = _mm256_add_ps(_mm256_mul_ps(xi_232, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_233, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)));
          const __m256 xi_240 = _mm256_add_ps(xi_232, xi_233);
          const __m256 xi_243 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_186, xi_194), xi_209), xi_227), xi_228), xi_229), xi_230);
          const __m256 forceTerm_0 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_25, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_26, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_27, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_25, _mm256_set_ps(xi_28, xi_28, xi_28, xi_28, xi_28, xi_28, xi_28, xi_28))), _mm256_mul_ps(xi_26, _mm256_set_ps(xi_28, xi_28, xi_28, xi_28, xi_28, xi_28, xi_28, xi_28))), _mm256_mul_ps(xi_27, _mm256_set_ps(xi_28, xi_28, xi_28, xi_28, xi_28, xi_28, xi_28, xi_28)));
          const __m256 forceTerm_1 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_31, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_29), xi_41);
          const __m256 forceTerm_2 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_29, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_31), xi_41);
          const __m256 forceTerm_3 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_42, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_44), xi_48);
          const __m256 forceTerm_4 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_44, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_42), xi_48);
          const __m256 forceTerm_5 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_51, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_49), xi_52);
          const __m256 forceTerm_6 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_49, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_51), xi_52);
          const __m256 forceTerm_7 = _mm256_add_ps(_mm256_add_ps(xi_57, xi_66), xi_76);
          const __m256 forceTerm_8 = _mm256_add_ps(_mm256_add_ps(xi_66, xi_77), xi_78);
          const __m256 forceTerm_9 = _mm256_add_ps(_mm256_add_ps(xi_57, xi_78), xi_80);
          const __m256 forceTerm_10 = _mm256_add_ps(_mm256_add_ps(xi_76, xi_77), xi_80);
          const __m256 forceTerm_11 = _mm256_add_ps(_mm256_add_ps(xi_83, xi_85), xi_90);
          const __m256 forceTerm_12 = _mm256_add_ps(_mm256_add_ps(xi_83, xi_91), xi_92);
          const __m256 forceTerm_13 = _mm256_add_ps(_mm256_add_ps(xi_101, xi_57), xi_94);
          const __m256 forceTerm_14 = _mm256_add_ps(_mm256_add_ps(xi_102, xi_77), xi_94);
          const __m256 forceTerm_15 = _mm256_add_ps(_mm256_add_ps(xi_103, xi_85), xi_92);
          const __m256 forceTerm_16 = _mm256_add_ps(_mm256_add_ps(xi_103, xi_90), xi_91);
          const __m256 forceTerm_17 = _mm256_add_ps(_mm256_add_ps(xi_102, xi_104), xi_57);
          const __m256 forceTerm_18 = _mm256_add_ps(_mm256_add_ps(xi_101, xi_104), xi_77);
          _mm256_store_ps(&_data_pdfs_20_30_10[ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_110, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_131, _mm256_set_ps(0.023809523809523808f, 0.023809523809523808f, 0.023809523809523808f, 0.023809523809523808f, 0.023809523809523808f, 0.023809523809523808f, 0.023809523809523808f, 0.023809523809523808f))), _mm256_mul_ps(xi_107, _mm256_set_ps(0.14285714285714285f, 0.14285714285714285f, 0.14285714285714285f, 0.14285714285714285f, 0.14285714285714285f, 0.14285714285714285f, 0.14285714285714285f, 0.14285714285714285f))), _mm256_mul_ps(xi_127, _mm256_set_ps(0.042857142857142858f, 0.042857142857142858f, 0.042857142857142858f, 0.042857142857142858f, 0.042857142857142858f, 0.042857142857142858f, 0.042857142857142858f, 0.042857142857142858f))), _mm256_mul_ps(xi_111, _mm256_set_ps(0.085714285714285715f, 0.085714285714285715f, 0.085714285714285715f, 0.085714285714285715f, 0.085714285714285715f, 0.085714285714285715f, 0.085714285714285715f, 0.085714285714285715f))), _mm256_mul_ps(xi_117, _mm256_set_ps(0.10000000000000001f, 0.10000000000000001f, 0.10000000000000001f, 0.10000000000000001f, 0.10000000000000001f, 0.10000000000000001f, 0.10000000000000001f, 0.10000000000000001f))), _mm256_mul_ps(xi_108, _mm256_set_ps(0.20000000000000001f, 0.20000000000000001f, 0.20000000000000001f, 0.20000000000000001f, 0.20000000000000001f, 0.20000000000000001f, 0.20000000000000001f, 0.20000000000000001f))), _mm256_mul_ps(xi_124, _mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f))), forceTerm_0), xi_251));
          _mm256_store_ps(&_data_pdfs_20_31_10[ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_135, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_140, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), forceTerm_1), xi_151), xi_164), xi_246));
          _mm256_store_ps(&_data_pdfs_20_32_10[ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_2, xi_135), xi_140), xi_164), xi_165), xi_247));
          _mm256_store_ps(&_data_pdfs_20_33_10[ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_3, xi_168), xi_170), xi_175), xi_176), xi_253));
          _mm256_store_ps(&_data_pdfs_20_34_10[ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_168, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_170, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), forceTerm_4), xi_176), xi_177), xi_252));
          _mm256_store_ps(&_data_pdfs_20_35_10[ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_179, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_181, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), forceTerm_5), xi_186), xi_187), xi_256));
          _mm256_store_ps(&_data_pdfs_20_36_10[ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_6, xi_179), xi_181), xi_187), xi_188), xi_244));
          _mm256_store_ps(&_data_pdfs_20_37_10[ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_7, xi_177), xi_195), xi_200), xi_265));
          _mm256_store_ps(&_data_pdfs_20_38_10[ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_8, xi_175), xi_195), xi_201), xi_260));
          _mm256_store_ps(&_data_pdfs_20_39_10[ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_9, xi_177), xi_201), xi_202), xi_264));
          _mm256_store_ps(&_data_pdfs_20_310_10[ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_10, xi_175), xi_200), xi_202), xi_263));
          _mm256_store_ps(&_data_pdfs_20_311_10[ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_11, xi_210), xi_219), xi_224), xi_248));
          _mm256_store_ps(&_data_pdfs_20_312_10[ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_12, xi_219), xi_225), xi_226), xi_258));
          _mm256_store_ps(&_data_pdfs_20_313_10[ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_13, xi_231), xi_234), xi_239), xi_250));
          _mm256_store_ps(&_data_pdfs_20_314_10[ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_14, xi_231), xi_240), xi_241), xi_257));
          _mm256_store_ps(&_data_pdfs_20_315_10[ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_15, xi_224), xi_226), xi_242), xi_261));
          _mm256_store_ps(&_data_pdfs_20_316_10[ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_16, xi_210), xi_225), xi_242), xi_259));
          _mm256_store_ps(&_data_pdfs_20_317_10[ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_17, xi_239), xi_240), xi_243), xi_254));
          _mm256_store_ps(&_data_pdfs_20_318_10[ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_18, xi_234), xi_241), xi_243), xi_249));
        }
        for (int64_t ctr_0 = (int64_t)((_size_force_0) / (8)) * (8); ctr_0 < _size_force_0; ctr_0 += 1) {
          const float xi_244 = _data_pdfs_20_36_10[ctr_0];
          const float xi_245 = _data_force_20_32_10[ctr_0];
          const float xi_246 = _data_pdfs_20_31_10[ctr_0];
          const float xi_247 = _data_pdfs_20_32_10[ctr_0];
          const float xi_248 = _data_pdfs_20_311_10[ctr_0];
          const float xi_249 = _data_pdfs_20_318_10[ctr_0];
          const float xi_250 = _data_pdfs_20_313_10[ctr_0];
          const float xi_251 = _data_pdfs_20_30_10[ctr_0];
          const float xi_252 = _data_pdfs_20_34_10[ctr_0];
          const float xi_253 = _data_pdfs_20_33_10[ctr_0];
          const float xi_254 = _data_pdfs_20_317_10[ctr_0];
          const float xi_255 = _data_force_20_30_10[ctr_0];
          const float xi_256 = _data_pdfs_20_35_10[ctr_0];
          const float xi_257 = _data_pdfs_20_314_10[ctr_0];
          const float xi_258 = _data_pdfs_20_312_10[ctr_0];
          const float xi_259 = _data_pdfs_20_316_10[ctr_0];
          const float xi_260 = _data_pdfs_20_38_10[ctr_0];
          const float xi_261 = _data_pdfs_20_315_10[ctr_0];
          const float xi_262 = _data_force_20_31_10[ctr_0];
          const float xi_263 = _data_pdfs_20_310_10[ctr_0];
          const float xi_264 = _data_pdfs_20_39_10[ctr_0];
          const float xi_265 = _data_pdfs_20_37_10[ctr_0];

          float random_3_0{};
          float random_3_1{};
          float random_3_2{};
          float random_3_3{};
          if (kT > 0.) {
            philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 3, seed, random_3_0, random_3_1, random_3_2, random_3_3);
          }

          float random_2_0{};
          float random_2_1{};
          float random_2_2{};
          float random_2_3{};
          if (kT > 0.) {
            philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 2, seed, random_2_0, random_2_1, random_2_2, random_2_3);
          }

          float random_1_0{};
          float random_1_1{};
          float random_1_2{};
          float random_1_3{};
          if (kT > 0.) {
            philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 1, seed, random_1_0, random_1_1, random_1_2, random_1_3);
          }

          float random_0_0{};
          float random_0_1{};
          float random_0_2{};
          float random_0_3{};
          if (kT > 0.) {
            philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 0, seed, random_0_0, random_0_1, random_0_2, random_0_3);
          }
          const float xi_2 = xi_249 + xi_257;
          const float xi_3 = xi_2 + xi_252;
          const float xi_4 = xi_246 + xi_248 + xi_261;
          const float xi_5 = xi_256 + xi_258;
          const float xi_6 = xi_244 + xi_254;
          const float xi_8 = xi_264 * -1.0f;
          const float xi_9 = xi_265 * -1.0f;
          const float xi_10 = xi_254 * -1.0f;
          const float xi_11 = xi_250 * -1.0f;
          const float xi_12 = xi_253 * -1.0f;
          const float xi_13 = xi_10 + xi_11 + xi_12;
          const float xi_14 = xi_247 * -1.0f;
          const float xi_15 = xi_263 * -1.0f;
          const float xi_16 = xi_14 + xi_15;
          const float xi_17 = xi_259 * -1.0f;
          const float xi_18 = xi_258 * -1.0f;
          const float xi_19 = xi_17 + xi_18;
          const float xi_20 = xi_249 * -1.0f;
          const float xi_21 = xi_10 + xi_20;
          const float xi_22 = xi_261 * -1.0f;
          const float xi_23 = xi_244 * -1.0f;
          const float xi_24 = xi_17 + xi_22 + xi_23 + xi_248;
          const float xi_29 = xi_262 * 0.16666666666666666f;
          const float xi_30 = xi_262 * 0.083333333333333329f;
          const float xi_42 = xi_255 * 0.16666666666666666f;
          const float xi_43 = xi_255 * 0.083333333333333329f;
          const float xi_49 = xi_245 * 0.16666666666666666f;
          const float xi_50 = xi_245 * 0.083333333333333329f;
          const float xi_67 = xi_262 * 0.25f;
          const float xi_72 = xi_262 * xi_71;
          const float xi_114 = xi_251 * -1.0f;
          const float xi_118 = xi_248 * -1.0f;
          const float xi_119 = xi_118 + xi_18;
          const float xi_120 = xi_260 * -1.0f + xi_8;
          const float xi_122 = xi_257 * -1.0f;
          const float xi_123 = xi_11 + xi_122 + xi_15 + xi_21;
          const float xi_125 = xi_248 * 2.0f + xi_258 * 2.0f + xi_259 * 2.0f + xi_261 * 2.0f;
          const float xi_126 = xi_125 + xi_252 * 5.0f + xi_253 * 5.0f;
          const float xi_128 = xi_250 * 2.0f;
          const float xi_129 = xi_257 * 2.0f;
          const float xi_130 = xi_249 * 2.0f + xi_254 * 2.0f;
          const float xi_132 = xi_118 + xi_258;
          const float xi_133 = xi_132 + xi_14 + xi_22 + xi_246 + xi_259;
          const float xi_135 = xi_133 * xi_134;
          const float xi_136 = random_2_3 - 0.5f;
          const float xi_141 = xi_265 * 2.0f;
          const float xi_142 = xi_263 * 2.0f;
          const float xi_143 = xi_260 * -2.0f + xi_264 * 2.0f;
          const float xi_144 = xi_14 + xi_141 * -1.0f + xi_142 + xi_143 + xi_19 + xi_4;
          const float xi_146 = xi_144 * xi_145;
          const float xi_147 = random_1_2 - 0.5f;
          const float xi_152 = random_0_1 - 0.5f;
          const float xi_166 = xi_122 + xi_250;
          const float xi_167 = xi_12 + xi_166 + xi_20 + xi_252 + xi_254;
          const float xi_168 = xi_134 * xi_167;
          const float xi_169 = random_2_1 - 0.5f;
          const float xi_171 = xi_13 + xi_141 + xi_142 * -1.0f + xi_143 + xi_3;
          const float xi_172 = xi_145 * xi_171;
          const float xi_173 = random_2_0 - 0.5f;
          const float xi_178 = xi_119 + xi_23 + xi_256 + xi_259 + xi_261;
          const float xi_179 = xi_134 * xi_178;
          const float xi_180 = random_2_2 - 0.5f;
          const float xi_182 = xi_128 * -1.0f + xi_129 * -1.0f + xi_130 + xi_24 + xi_5;
          const float xi_183 = xi_145 * xi_182;
          const float xi_184 = random_1_3 - 0.5f;
          const float xi_212 = xi_182 * xi_211;
          const float xi_214 = xi_178 * xi_213;
          const float xi_220 = xi_144 * xi_211;
          const float xi_221 = xi_133 * xi_213;
          const float xi_235 = xi_167 * xi_213;
          const float xi_236 = xi_171 * xi_211;
          const float xi_31 = rr_0 * xi_30;
          const float xi_44 = rr_0 * xi_43;
          const float xi_51 = rr_0 * xi_50;
          const float xi_54 = xi_255 * xi_53;
          const float xi_59 = xi_262 * xi_53;
          const float xi_81 = xi_245 * xi_53;
          const float vel0Term = xi_260 + xi_263 + xi_3;
          const float vel1Term = xi_265 + xi_4;
          const float vel2Term = xi_250 + xi_5;
          const float rho = vel0Term + vel1Term + vel2Term + xi_247 + xi_251 + xi_253 + xi_259 + xi_264 + xi_6;
          const float xi_105 = kT * rho;
          const float xi_106 = powf(xi_105 * (-1.0f * (omega_even * -1.0f + 1.0f) * (omega_even * -1.0f + 1.0f) + 1.0f), 0.5f);
          const float xi_107 = xi_106 * (random_3_0 - 0.5f) * 3.7416573867739413f;
          const float xi_108 = xi_106 * (random_3_2 - 0.5f) * 5.4772255750516612f;
          const float xi_110 = xi_109 * (random_1_1 - 0.5f) * powf(xi_105 * (-1.0f * (omega_bulk * -1.0f + 1.0f) * (omega_bulk * -1.0f + 1.0f) + 1.0f), 0.5f);
          const float xi_111 = xi_106 * (random_3_1 - 0.5f) * 8.3666002653407556f;
          const float xi_137 = powf(xi_105 * (-1.0f * (omega_odd * -1.0f + 1.0f) * (omega_odd * -1.0f + 1.0f) + 1.0f), 0.5f);
          const float xi_138 = xi_137 * 1.4142135623730951f;
          const float xi_139 = xi_138 * 0.5f;
          const float xi_140 = xi_136 * xi_139;
          const float xi_148 = xi_109 * xi_137;
          const float xi_149 = xi_148 * 0.16666666666666666f;
          const float xi_150 = xi_147 * xi_149;
          const float xi_151 = xi_146 * -1.0f + xi_150 * -1.0f;
          const float xi_153 = powf(xi_105 * (-1.0f * (omega_shear * -1.0f + 1.0f) * (omega_shear * -1.0f + 1.0f) + 1.0f), 0.5f);
          const float xi_154 = xi_153 * 0.5f;
          const float xi_155 = xi_152 * xi_154;
          const float xi_161 = xi_153 * (random_0_0 - 0.5f) * 1.7320508075688772f;
          const float xi_165 = xi_146 + xi_150;
          const float xi_170 = xi_139 * xi_169;
          const float xi_174 = xi_149 * xi_173;
          const float xi_175 = xi_172 + xi_174;
          const float xi_177 = xi_172 * -1.0f + xi_174 * -1.0f;
          const float xi_181 = xi_139 * xi_180;
          const float xi_185 = xi_149 * xi_184;
          const float xi_186 = xi_183 * -1.0f + xi_185 * -1.0f;
          const float xi_188 = xi_183 + xi_185;
          const float xi_189 = xi_152 * xi_153 * 0.25f;
          const float xi_192 = xi_107 * 0.083333333333333329f;
          const float xi_196 = xi_154 * (random_0_2 - 0.5f);
          const float xi_203 = xi_154 * (random_1_0 - 0.5f);
          const float xi_207 = xi_111 * -0.014285714285714285f;
          const float xi_208 = xi_108 * 0.050000000000000003f;
          const float xi_215 = xi_148 * 0.083333333333333329f;
          const float xi_216 = xi_184 * xi_215;
          const float xi_217 = xi_138 * 0.25f;
          const float xi_218 = xi_180 * xi_217;
          const float xi_219 = xi_212 * -1.0f + xi_214 + xi_216 * -1.0f + xi_218;
          const float xi_222 = xi_147 * xi_215;
          const float xi_223 = xi_136 * xi_217;
          const float xi_224 = xi_220 * -1.0f + xi_221 + xi_222 * -1.0f + xi_223;
          const float xi_225 = xi_220 + xi_221 * -1.0f + xi_222 + xi_223 * -1.0f;
          const float xi_227 = xi_189 * -1.0f;
          const float xi_230 = xi_111 * 0.035714285714285712f;
          const float xi_232 = xi_154 * (random_0_3 - 0.5f);
          const float xi_237 = xi_169 * xi_217;
          const float xi_238 = xi_173 * xi_215;
          const float xi_239 = xi_235 * -1.0f + xi_236 + xi_237 * -1.0f + xi_238;
          const float xi_241 = xi_235 + xi_236 * -1.0f + xi_237 + xi_238 * -1.0f;
          const float xi_242 = xi_212 + xi_214 * -1.0f + xi_216 + xi_218 * -1.0f;
          const float xi_0 = ((1.0f) / (rho));
          const float xi_7 = xi_0 * 0.5f;
          const float u_0 = xi_0 * (vel0Term + xi_13 + xi_8 + xi_9) + xi_255 * xi_7;
          const float xi_25 = u_0 * xi_255;
          const float xi_37 = xi_25 * 0.16666666666666666f;
          const float xi_38 = xi_25 * 0.083333333333333329f;
          const float xi_39 = omega_shear * xi_38;
          const float xi_40 = xi_37 * -1.0f + xi_39;
          const float xi_56 = xi_25 * xi_55 * -1.0f + xi_37;
          const float xi_57 = xi_43 * -1.0f + xi_54 + xi_56;
          const float xi_61 = xi_25 * xi_60 * -1.0f;
          const float xi_68 = u_0 * xi_67;
          const float xi_73 = u_0 * xi_72;
          const float xi_77 = xi_43 + xi_54 * -1.0f + xi_56;
          const float xi_84 = xi_38 * -1.0f;
          const float xi_95 = u_0 * xi_245;
          const float xi_96 = xi_95 * 0.25f;
          const float xi_99 = xi_71 * xi_95;
          const float xi_113 = rho * u_0 * u_0;
          const float u_1 = xi_0 * (vel1Term + xi_16 + xi_19 + xi_260 + xi_8) + xi_262 * xi_7;
          const float xi_26 = u_1 * xi_262;
          const float xi_32 = xi_26 * 0.16666666666666666f;
          const float xi_45 = xi_26 * 0.083333333333333329f;
          const float xi_46 = omega_shear * xi_45;
          const float xi_47 = xi_32 * -1.0f + xi_46;
          const float xi_62 = xi_26 * xi_60 * -1.0f;
          const float xi_69 = u_1 * 0.25f;
          const float xi_70 = xi_255 * xi_69;
          const float xi_74 = u_1 * xi_71;
          const float xi_75 = xi_255 * xi_74;
          const float xi_76 = xi_68 * -1.0f + xi_70 * -1.0f + xi_73 + xi_75;
          const float xi_78 = xi_68 + xi_70 + xi_73 * -1.0f + xi_75 * -1.0f;
          const float xi_86 = xi_245 * xi_69;
          const float xi_88 = xi_245 * xi_74;
          const float xi_93 = xi_45 * -1.0f;
          const float xi_112 = rho * u_1 * u_1;
          const float xi_121 = xi_112 + xi_120 + xi_9;
          const float xi_197 = rho * u_1;
          const float xi_199 = xi_198 * (u_0 * xi_197 + xi_120 + xi_263 + xi_265);
          const float xi_200 = xi_196 * -1.0f + xi_199 * -1.0f;
          const float xi_201 = xi_196 + xi_199;
          const float u_2 = xi_0 * (vel2Term + xi_21 + xi_24 + xi_257) + xi_245 * xi_7;
          const float xi_27 = u_2 * xi_245;
          const float xi_33 = xi_27 * 0.16666666666666666f;
          const float xi_34 = xi_27 * 0.083333333333333329f;
          const float xi_35 = omega_shear * xi_34;
          const float xi_36 = xi_33 * -1.0f + xi_35;
          const float xi_41 = omega_shear * xi_32 * -1.0f + xi_26 * 0.33333333333333331f + xi_36 + xi_40;
          const float xi_48 = omega_shear * xi_37 * -1.0f + xi_25 * 0.33333333333333331f + xi_36 + xi_47;
          const float xi_52 = omega_shear * xi_33 * -1.0f + xi_27 * 0.33333333333333331f + xi_40 + xi_47;
          const float xi_58 = xi_34 * -1.0f;
          const float xi_63 = xi_27 * xi_60 * -1.0f;
          const float xi_64 = xi_26 * xi_55 * -1.0f + xi_32 + xi_61 + xi_62 + xi_63;
          const float xi_65 = xi_30 + xi_59 * -1.0f + xi_64;
          const float xi_66 = xi_35 + xi_58 + xi_65;
          const float xi_79 = xi_30 * -1.0f + xi_59 + xi_64;
          const float xi_80 = xi_35 + xi_58 + xi_79;
          const float xi_82 = xi_27 * xi_55 * -1.0f + xi_33;
          const float xi_83 = xi_50 + xi_81 * -1.0f + xi_82;
          const float xi_85 = xi_39 + xi_65 + xi_84;
          const float xi_87 = u_2 * xi_67;
          const float xi_89 = u_2 * xi_72;
          const float xi_90 = xi_86 + xi_87 + xi_88 * -1.0f + xi_89 * -1.0f;
          const float xi_91 = xi_39 + xi_79 + xi_84;
          const float xi_92 = xi_86 * -1.0f + xi_87 * -1.0f + xi_88 + xi_89;
          const float xi_94 = xi_46 + xi_61 + xi_62 + xi_63 + xi_83 + xi_93;
          const float xi_97 = u_2 * xi_255;
          const float xi_98 = xi_97 * 0.25f;
          const float xi_100 = xi_71 * xi_97;
          const float xi_101 = xi_100 + xi_96 * -1.0f + xi_98 * -1.0f + xi_99;
          const float xi_102 = xi_100 * -1.0f + xi_96 + xi_98 + xi_99 * -1.0f;
          const float xi_103 = xi_50 * -1.0f + xi_81 + xi_82;
          const float xi_104 = xi_103 + xi_46 + xi_61 + xi_62 + xi_63 + xi_93;
          const float xi_115 = rho * u_2 * u_2;
          const float xi_116 = xi_114 + xi_115 * 0.66666666666666663f + xi_244 * 3.0f + xi_256 * 3.0f;
          const float xi_117 = omega_even * (xi_112 * 0.66666666666666663f + xi_113 * 1.6666666666666667f + xi_116 + xi_246 * 3.0f + xi_247 * 3.0f + xi_248 * -3.0f + xi_258 * -3.0f + xi_259 * -3.0f + xi_261 * -3.0f);
          const float xi_124 = omega_bulk * (xi_113 + xi_115 + xi_119 + xi_121 + xi_123 + xi_17 + xi_22 + xi_251);
          const float xi_127 = omega_even * (xi_112 * 2.3333333333333335f + xi_116 + xi_126 + xi_246 * -2.0f + xi_247 * -2.0f + xi_249 * -5.0f + xi_250 * -5.0f + xi_254 * -5.0f + xi_257 * -5.0f);
          const float xi_131 = omega_even * (xi_114 + xi_115 * 3.0f + xi_126 + xi_128 + xi_129 + xi_130 + xi_244 * -4.0f + xi_246 * 5.0f + xi_247 * 5.0f + xi_256 * -4.0f + xi_260 * -7.0f + xi_263 * -7.0f + xi_264 * -7.0f + xi_265 * -7.0f);
          const float xi_156 = xi_115 * -1.0f + xi_256;
          const float xi_157 = omega_shear * (xi_121 + xi_156 + xi_16 + xi_2 + xi_246 * -1.0f + xi_250 + xi_6);
          const float xi_158 = xi_157 * 0.125f;
          const float xi_159 = xi_107 * -0.11904761904761904f + xi_131 * -0.01984126984126984f;
          const float xi_160 = omega_shear * (xi_112 * -1.0f + xi_113 * 2.0f + xi_120 + xi_123 + xi_125 + xi_156 + xi_244 + xi_246 + xi_247 + xi_252 * -2.0f + xi_253 * -2.0f + xi_9);
          const float xi_162 = xi_160 * -0.041666666666666664f + xi_161 * -0.16666666666666666f;
          const float xi_163 = xi_108 * -0.10000000000000001f + xi_117 * -0.050000000000000003f + xi_162;
          const float xi_164 = xi_111 * 0.028571428571428571f + xi_127 * 0.014285714285714285f + xi_155 + xi_158 + xi_159 + xi_163;
          const float xi_176 = xi_111 * -0.071428571428571425f + xi_127 * -0.035714285714285712f + xi_159 + xi_160 * 0.083333333333333329f + xi_161 * 0.33333333333333331f;
          const float xi_187 = xi_107 * 0.095238095238095233f + xi_111 * -0.042857142857142858f + xi_127 * -0.021428571428571429f + xi_131 * 0.015873015873015872f + xi_155 * -1.0f + xi_158 * -1.0f + xi_163;
          const float xi_190 = xi_157 * 0.0625f;
          const float xi_191 = xi_131 * 0.013888888888888888f;
          const float xi_193 = xi_110 * 0.083333333333333329f + xi_124 * 0.041666666666666664f;
          const float xi_194 = xi_160 * 0.020833333333333332f + xi_161 * 0.083333333333333329f + xi_193;
          const float xi_195 = xi_165 + xi_189 + xi_190 + xi_191 + xi_192 + xi_194;
          const float xi_202 = xi_151 + xi_189 + xi_190 + xi_191 + xi_192 + xi_194;
          const float xi_204 = xi_127 * -0.0071428571428571426f;
          const float xi_205 = xi_198 * (u_2 * xi_197 + xi_132 + xi_17 + xi_261);
          const float xi_206 = xi_117 * 0.025000000000000001f;
          const float xi_209 = xi_107 * -0.023809523809523808f + xi_131 * -0.003968253968253968f;
          const float xi_210 = xi_162 + xi_193 + xi_203 + xi_204 + xi_205 + xi_206 + xi_207 + xi_208 + xi_209;
          const float xi_226 = xi_162 + xi_193 + xi_203 * -1.0f + xi_204 + xi_205 * -1.0f + xi_206 + xi_207 + xi_208 + xi_209;
          const float xi_228 = xi_190 * -1.0f;
          const float xi_229 = xi_127 * 0.017857142857142856f;
          const float xi_231 = xi_188 + xi_194 + xi_209 + xi_227 + xi_228 + xi_229 + xi_230;
          const float xi_233 = xi_198 * (rho * u_0 * u_2 + xi_10 + xi_166 + xi_249);
          const float xi_234 = xi_232 * -1.0f + xi_233 * -1.0f;
          const float xi_240 = xi_232 + xi_233;
          const float xi_243 = xi_186 + xi_194 + xi_209 + xi_227 + xi_228 + xi_229 + xi_230;
          const float forceTerm_0 = xi_25 * xi_28 + xi_25 * -1.0f + xi_26 * xi_28 + xi_26 * -1.0f + xi_27 * xi_28 + xi_27 * -1.0f;
          const float forceTerm_1 = xi_29 + xi_31 * -1.0f + xi_41;
          const float forceTerm_2 = xi_29 * -1.0f + xi_31 + xi_41;
          const float forceTerm_3 = xi_42 * -1.0f + xi_44 + xi_48;
          const float forceTerm_4 = xi_42 + xi_44 * -1.0f + xi_48;
          const float forceTerm_5 = xi_49 + xi_51 * -1.0f + xi_52;
          const float forceTerm_6 = xi_49 * -1.0f + xi_51 + xi_52;
          const float forceTerm_7 = xi_57 + xi_66 + xi_76;
          const float forceTerm_8 = xi_66 + xi_77 + xi_78;
          const float forceTerm_9 = xi_57 + xi_78 + xi_80;
          const float forceTerm_10 = xi_76 + xi_77 + xi_80;
          const float forceTerm_11 = xi_83 + xi_85 + xi_90;
          const float forceTerm_12 = xi_83 + xi_91 + xi_92;
          const float forceTerm_13 = xi_101 + xi_57 + xi_94;
          const float forceTerm_14 = xi_102 + xi_77 + xi_94;
          const float forceTerm_15 = xi_103 + xi_85 + xi_92;
          const float forceTerm_16 = xi_103 + xi_90 + xi_91;
          const float forceTerm_17 = xi_102 + xi_104 + xi_57;
          const float forceTerm_18 = xi_101 + xi_104 + xi_77;
          _data_pdfs_20_30_10[ctr_0] = forceTerm_0 + xi_107 * 0.14285714285714285f + xi_108 * 0.20000000000000001f + xi_110 * -1.0f + xi_111 * 0.085714285714285715f + xi_117 * 0.10000000000000001f + xi_124 * -0.5f + xi_127 * 0.042857142857142858f + xi_131 * 0.023809523809523808f + xi_251;
          _data_pdfs_20_31_10[ctr_0] = forceTerm_1 + xi_135 * -1.0f + xi_140 * -1.0f + xi_151 + xi_164 + xi_246;
          _data_pdfs_20_32_10[ctr_0] = forceTerm_2 + xi_135 + xi_140 + xi_164 + xi_165 + xi_247;
          _data_pdfs_20_33_10[ctr_0] = forceTerm_3 + xi_168 + xi_170 + xi_175 + xi_176 + xi_253;
          _data_pdfs_20_34_10[ctr_0] = forceTerm_4 + xi_168 * -1.0f + xi_170 * -1.0f + xi_176 + xi_177 + xi_252;
          _data_pdfs_20_35_10[ctr_0] = forceTerm_5 + xi_179 * -1.0f + xi_181 * -1.0f + xi_186 + xi_187 + xi_256;
          _data_pdfs_20_36_10[ctr_0] = forceTerm_6 + xi_179 + xi_181 + xi_187 + xi_188 + xi_244;
          _data_pdfs_20_37_10[ctr_0] = forceTerm_7 + xi_177 + xi_195 + xi_200 + xi_265;
          _data_pdfs_20_38_10[ctr_0] = forceTerm_8 + xi_175 + xi_195 + xi_201 + xi_260;
          _data_pdfs_20_39_10[ctr_0] = forceTerm_9 + xi_177 + xi_201 + xi_202 + xi_264;
          _data_pdfs_20_310_10[ctr_0] = forceTerm_10 + xi_175 + xi_200 + xi_202 + xi_263;
          _data_pdfs_20_311_10[ctr_0] = forceTerm_11 + xi_210 + xi_219 + xi_224 + xi_248;
          _data_pdfs_20_312_10[ctr_0] = forceTerm_12 + xi_219 + xi_225 + xi_226 + xi_258;
          _data_pdfs_20_313_10[ctr_0] = forceTerm_13 + xi_231 + xi_234 + xi_239 + xi_250;
          _data_pdfs_20_314_10[ctr_0] = forceTerm_14 + xi_231 + xi_240 + xi_241 + xi_257;
          _data_pdfs_20_315_10[ctr_0] = forceTerm_15 + xi_224 + xi_226 + xi_242 + xi_261;
          _data_pdfs_20_316_10[ctr_0] = forceTerm_16 + xi_210 + xi_225 + xi_242 + xi_259;
          _data_pdfs_20_317_10[ctr_0] = forceTerm_17 + xi_239 + xi_240 + xi_243 + xi_254;
          _data_pdfs_20_318_10[ctr_0] = forceTerm_18 + xi_234 + xi_241 + xi_243 + xi_249;
        }
      }
    }
  }
}
} // namespace internal_48c9ee502281a70505dce0378c55abd5

void CollideSweepSinglePrecisionThermalizedAVX::run(IBlock *block) {
  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);
  auto force = block->getData<field::GhostLayerField<float, 3>>(forceID);

  auto &time_step = this->time_step_;
  auto &kT = this->kT_;
  auto &omega_odd = this->omega_odd_;
  auto &seed = this->seed_;
  auto &omega_bulk = this->omega_bulk_;
  auto block_offset_0 = this->block_offset_0_;
  auto &omega_shear = this->omega_shear_;
  auto &omega_even = this->omega_even_;
  auto block_offset_2 = this->block_offset_2_;
  auto block_offset_1 = this->block_offset_1_;
  block_offset_generator(block, block_offset_0, block_offset_1, block_offset_2);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force->nrOfGhostLayers()));
  float *RESTRICT const _data_force = force->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()));
  float *RESTRICT _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
  WALBERLA_ASSERT_EQUAL((uintptr_t)pdfs->dataAt(0, 0, 0, 0) % 32, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(), int64_t(cell_idx_c(force->xSize()) + 0));
  const int64_t _size_force_0 = int64_t(cell_idx_c(force->xSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(), int64_t(cell_idx_c(force->ySize()) + 0));
  const int64_t _size_force_1 = int64_t(cell_idx_c(force->ySize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(), int64_t(cell_idx_c(force->zSize()) + 0));
  const int64_t _size_force_2 = int64_t(cell_idx_c(force->zSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0);
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  internal_48c9ee502281a70505dce0378c55abd5::collidesweepsingleprecisionthermalizedavx_collidesweepsingleprecisionthermalizedavx(_data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, block_offset_0, block_offset_1, block_offset_2, kT, omega_bulk, omega_even, omega_odd, omega_shear, seed, time_step);
}

void CollideSweepSinglePrecisionThermalizedAVX::runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks, const CellInterval &globalCellInterval, cell_idx_t ghostLayers, IBlock *block) {
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
  auto &kT = this->kT_;
  auto &omega_odd = this->omega_odd_;
  auto &seed = this->seed_;
  auto &omega_bulk = this->omega_bulk_;
  auto block_offset_0 = this->block_offset_0_;
  auto &omega_shear = this->omega_shear_;
  auto &omega_even = this->omega_even_;
  auto block_offset_2 = this->block_offset_2_;
  auto block_offset_1 = this->block_offset_1_;
  block_offset_generator(block, block_offset_0, block_offset_1, block_offset_2);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(force->nrOfGhostLayers()));
  float *RESTRICT const _data_force = force->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
  float *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
  WALBERLA_ASSERT_EQUAL((uintptr_t)pdfs->dataAt(0, 0, 0, 0) % 32, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
  const int64_t _size_force_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
  const int64_t _size_force_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
  const int64_t _size_force_2 = int64_t(cell_idx_c(ci.zSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0);
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  internal_48c9ee502281a70505dce0378c55abd5::collidesweepsingleprecisionthermalizedavx_collidesweepsingleprecisionthermalizedavx(_data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, block_offset_0, block_offset_1, block_offset_2, kT, omega_bulk, omega_even, omega_odd, omega_shear, seed, time_step);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif