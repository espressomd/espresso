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
//! \\file CollideSweepDoublePrecisionThermalizedAVX.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

// kernel generated with pystencils v1.2, lbmpy v1.2, lbmpy_walberla/pystencils_walberla from waLBerla commit 4d10e7f2358fc4a4f7e99195d0f67f0b759ecb6f

#include <cmath>

#include "CollideSweepDoublePrecisionThermalizedAVX.h"
#include "core/DataTypes.h"
#include "core/Macros.h"

#include <immintrin.h>

#include "philox_rand.h"

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

namespace internal_25bc51f30ec2c20f3ee9796f7dcb65c6 {
static FUNC_PREFIX void collidesweepdoubleprecisionthermalizedavx_collidesweepdoubleprecisionthermalizedavx(double *RESTRICT const _data_force, double *RESTRICT _data_pdfs, int64_t const _size_force_0, int64_t const _size_force_1, int64_t const _size_force_2, int64_t const _stride_force_1, int64_t const _stride_force_2, int64_t const _stride_force_3, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, uint32_t block_offset_0, uint32_t block_offset_1, uint32_t block_offset_2, double kT, double omega_bulk, double omega_even, double omega_odd, double omega_shear, uint32_t seed, uint32_t time_step) {
  const double xi_28 = omega_bulk * 0.5;
  const double xi_55 = omega_shear * 0.041666666666666664;
  const double xi_60 = omega_bulk * 0.041666666666666664;
  const double xi_71 = omega_shear * 0.125;
  const double xi_109 = 2.4494897427831779;
  const double xi_134 = omega_odd * 0.25;
  const double xi_145 = omega_odd * 0.083333333333333329;
  const double xi_198 = omega_shear * 0.25;
  const double xi_211 = omega_odd * 0.041666666666666664;
  const double xi_213 = omega_odd * 0.125;
  const double rr_0 = 0.0;
  const double xi_53 = rr_0 * 0.041666666666666664;
  for (int64_t ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_20_34 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_36 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_315 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_310 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_312 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_318 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_39 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_31 = _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_37 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    double *RESTRICT _data_force_20_31 = _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    double *RESTRICT _data_pdfs_20_316 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_313 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_38 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_33 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    double *RESTRICT _data_force_20_32 = _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    double *RESTRICT _data_pdfs_20_314 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3;
    double *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    double *RESTRICT _data_pdfs_20_317 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_311 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_32 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_35 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_20_34_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      double *RESTRICT _data_pdfs_20_36_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_36;
      double *RESTRICT _data_pdfs_20_315_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_315;
      double *RESTRICT _data_pdfs_20_310_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_310;
      double *RESTRICT _data_pdfs_20_312_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_312;
      double *RESTRICT _data_pdfs_20_318_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_318;
      double *RESTRICT _data_pdfs_20_39_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_39;
      double *RESTRICT _data_pdfs_20_31_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_31;
      double *RESTRICT _data_pdfs_20_37_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_37;
      double *RESTRICT _data_pdfs_20_30_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      double *RESTRICT _data_force_20_31_10 = _stride_force_1 * ctr_1 + _data_force_20_31;
      double *RESTRICT _data_pdfs_20_316_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_316;
      double *RESTRICT _data_pdfs_20_313_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_313;
      double *RESTRICT _data_pdfs_20_38_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_38;
      double *RESTRICT _data_pdfs_20_33_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      double *RESTRICT _data_force_20_32_10 = _stride_force_1 * ctr_1 + _data_force_20_32;
      double *RESTRICT _data_pdfs_20_314_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_314;
      double *RESTRICT _data_force_20_30_10 = _stride_force_1 * ctr_1 + _data_force_20_30;
      double *RESTRICT _data_pdfs_20_317_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_317;
      double *RESTRICT _data_pdfs_20_311_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_311;
      double *RESTRICT _data_pdfs_20_32_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_32;
      double *RESTRICT _data_pdfs_20_35_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_35;
      {
        for (int64_t ctr_0 = 0; ctr_0 < (int64_t)((_size_force_0) / (4)) * (4); ctr_0 += 4) {
          const __m256d xi_244 = _mm256_load_pd(&_data_pdfs_20_34_10[ctr_0]);
          const __m256d xi_245 = _mm256_load_pd(&_data_pdfs_20_36_10[ctr_0]);
          const __m256d xi_246 = _mm256_load_pd(&_data_pdfs_20_315_10[ctr_0]);
          const __m256d xi_247 = _mm256_load_pd(&_data_pdfs_20_310_10[ctr_0]);
          const __m256d xi_248 = _mm256_load_pd(&_data_pdfs_20_312_10[ctr_0]);
          const __m256d xi_249 = _mm256_load_pd(&_data_pdfs_20_318_10[ctr_0]);
          const __m256d xi_250 = _mm256_load_pd(&_data_pdfs_20_39_10[ctr_0]);
          const __m256d xi_251 = _mm256_load_pd(&_data_pdfs_20_31_10[ctr_0]);
          const __m256d xi_252 = _mm256_load_pd(&_data_pdfs_20_37_10[ctr_0]);
          const __m256d xi_253 = _mm256_load_pd(&_data_pdfs_20_30_10[ctr_0]);
          const __m256d xi_254 = _mm256_load_pd(&_data_force_20_31_10[ctr_0]);
          const __m256d xi_255 = _mm256_load_pd(&_data_pdfs_20_316_10[ctr_0]);
          const __m256d xi_256 = _mm256_load_pd(&_data_pdfs_20_313_10[ctr_0]);
          const __m256d xi_257 = _mm256_load_pd(&_data_pdfs_20_38_10[ctr_0]);
          const __m256d xi_258 = _mm256_load_pd(&_data_pdfs_20_33_10[ctr_0]);
          const __m256d xi_259 = _mm256_load_pd(&_data_force_20_32_10[ctr_0]);
          const __m256d xi_260 = _mm256_load_pd(&_data_pdfs_20_314_10[ctr_0]);
          const __m256d xi_261 = _mm256_load_pd(&_data_force_20_30_10[ctr_0]);
          const __m256d xi_262 = _mm256_load_pd(&_data_pdfs_20_317_10[ctr_0]);
          const __m256d xi_263 = _mm256_load_pd(&_data_pdfs_20_311_10[ctr_0]);
          const __m256d xi_264 = _mm256_load_pd(&_data_pdfs_20_32_10[ctr_0]);
          const __m256d xi_265 = _mm256_load_pd(&_data_pdfs_20_35_10[ctr_0]);

          __m256d random_7_0{};
          __m256d random_7_1{};
          if (kT > 0.) {
            philox_double2(time_step, _mm256_add_epi32(_mm256_add_epi32(_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0), _mm256_set_epi32(ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0)), _mm256_set_epi32(((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)))), block_offset_1 + ctr_1, block_offset_2 + ctr_2, 7, seed, random_7_0, random_7_1);
          }

          __m256d random_6_0{};
          __m256d random_6_1{};
          if (kT > 0.) {
            philox_double2(time_step, _mm256_add_epi32(_mm256_add_epi32(_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0), _mm256_set_epi32(ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0)), _mm256_set_epi32(((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)))), block_offset_1 + ctr_1, block_offset_2 + ctr_2, 6, seed, random_6_0, random_6_1);
          }

          __m256d random_5_0{};
          __m256d random_5_1{};
          if (kT > 0.) {
            philox_double2(time_step, _mm256_add_epi32(_mm256_add_epi32(_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0), _mm256_set_epi32(ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0)), _mm256_set_epi32(((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)))), block_offset_1 + ctr_1, block_offset_2 + ctr_2, 5, seed, random_5_0, random_5_1);
          }

          __m256d random_4_0{};
          __m256d random_4_1{};
          if (kT > 0.) {
            philox_double2(time_step, _mm256_add_epi32(_mm256_add_epi32(_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0), _mm256_set_epi32(ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0)), _mm256_set_epi32(((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)))), block_offset_1 + ctr_1, block_offset_2 + ctr_2, 4, seed, random_4_0, random_4_1);
          }

          __m256d random_3_0{};
          __m256d random_3_1{};
          if (kT > 0.) {
            philox_double2(time_step, _mm256_add_epi32(_mm256_add_epi32(_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0), _mm256_set_epi32(ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0)), _mm256_set_epi32(((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)))), block_offset_1 + ctr_1, block_offset_2 + ctr_2, 3, seed, random_3_0, random_3_1);
          }

          __m256d random_2_0{};
          __m256d random_2_1{};
          if (kT > 0.) {
            philox_double2(time_step, _mm256_add_epi32(_mm256_add_epi32(_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0), _mm256_set_epi32(ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0)), _mm256_set_epi32(((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)))), block_offset_1 + ctr_1, block_offset_2 + ctr_2, 2, seed, random_2_0, random_2_1);
          }

          __m256d random_1_0{};
          __m256d random_1_1{};
          if (kT > 0.) {
            philox_double2(time_step, _mm256_add_epi32(_mm256_add_epi32(_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0), _mm256_set_epi32(ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0)), _mm256_set_epi32(((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)))), block_offset_1 + ctr_1, block_offset_2 + ctr_2, 1, seed, random_1_0, random_1_1);
          }

          __m256d random_0_0{};
          __m256d random_0_1{};
          if (kT > 0.) {
            philox_double2(time_step, _mm256_add_epi32(_mm256_add_epi32(_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0), _mm256_set_epi32(ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0, ctr_0)), _mm256_set_epi32(((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)), ((int64_t)(block_offset_0)))), block_offset_1 + ctr_1, block_offset_2 + ctr_2, 0, seed, random_0_0, random_0_1);
          }
          const __m256d xi_2 = _mm256_add_pd(xi_249, xi_260);
          const __m256d xi_3 = _mm256_add_pd(xi_2, xi_244);
          const __m256d xi_4 = _mm256_add_pd(_mm256_add_pd(xi_246, xi_251), xi_263);
          const __m256d xi_5 = _mm256_add_pd(xi_248, xi_265);
          const __m256d xi_6 = _mm256_add_pd(xi_245, xi_262);
          const __m256d xi_8 = _mm256_mul_pd(xi_250, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_9 = _mm256_mul_pd(xi_252, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_10 = _mm256_mul_pd(xi_262, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_11 = _mm256_mul_pd(xi_256, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_12 = _mm256_mul_pd(xi_258, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_13 = _mm256_add_pd(_mm256_add_pd(xi_10, xi_11), xi_12);
          const __m256d xi_14 = _mm256_mul_pd(xi_264, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_15 = _mm256_mul_pd(xi_247, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_16 = _mm256_add_pd(xi_14, xi_15);
          const __m256d xi_17 = _mm256_mul_pd(xi_255, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_18 = _mm256_mul_pd(xi_248, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_19 = _mm256_add_pd(xi_17, xi_18);
          const __m256d xi_20 = _mm256_mul_pd(xi_249, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_21 = _mm256_add_pd(xi_10, xi_20);
          const __m256d xi_22 = _mm256_mul_pd(xi_246, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_23 = _mm256_mul_pd(xi_245, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_24 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_17, xi_22), xi_23), xi_263);
          const __m256d xi_29 = _mm256_mul_pd(xi_254, _mm256_set_pd(0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666666));
          const __m256d xi_30 = _mm256_mul_pd(xi_254, _mm256_set_pd(0.083333333333333329, 0.083333333333333329, 0.083333333333333329, 0.083333333333333329));
          const __m256d xi_42 = _mm256_mul_pd(xi_261, _mm256_set_pd(0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666666));
          const __m256d xi_43 = _mm256_mul_pd(xi_261, _mm256_set_pd(0.083333333333333329, 0.083333333333333329, 0.083333333333333329, 0.083333333333333329));
          const __m256d xi_49 = _mm256_mul_pd(xi_259, _mm256_set_pd(0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666666));
          const __m256d xi_50 = _mm256_mul_pd(xi_259, _mm256_set_pd(0.083333333333333329, 0.083333333333333329, 0.083333333333333329, 0.083333333333333329));
          const __m256d xi_67 = _mm256_mul_pd(xi_254, _mm256_set_pd(0.25, 0.25, 0.25, 0.25));
          const __m256d xi_72 = _mm256_mul_pd(xi_254, _mm256_set_pd(xi_71, xi_71, xi_71, xi_71));
          const __m256d xi_114 = _mm256_mul_pd(xi_253, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_118 = _mm256_mul_pd(xi_263, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_119 = _mm256_add_pd(xi_118, xi_18);
          const __m256d xi_120 = _mm256_add_pd(_mm256_mul_pd(xi_257, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_8);
          const __m256d xi_122 = _mm256_mul_pd(xi_260, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_123 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_11, xi_122), xi_15), xi_21);
          const __m256d xi_125 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_246, _mm256_set_pd(2.0, 2.0, 2.0, 2.0)), _mm256_mul_pd(xi_248, _mm256_set_pd(2.0, 2.0, 2.0, 2.0))), _mm256_mul_pd(xi_255, _mm256_set_pd(2.0, 2.0, 2.0, 2.0))), _mm256_mul_pd(xi_263, _mm256_set_pd(2.0, 2.0, 2.0, 2.0)));
          const __m256d xi_126 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_244, _mm256_set_pd(5.0, 5.0, 5.0, 5.0)), _mm256_mul_pd(xi_258, _mm256_set_pd(5.0, 5.0, 5.0, 5.0))), xi_125);
          const __m256d xi_128 = _mm256_mul_pd(xi_256, _mm256_set_pd(2.0, 2.0, 2.0, 2.0));
          const __m256d xi_129 = _mm256_mul_pd(xi_260, _mm256_set_pd(2.0, 2.0, 2.0, 2.0));
          const __m256d xi_130 = _mm256_add_pd(_mm256_mul_pd(xi_249, _mm256_set_pd(2.0, 2.0, 2.0, 2.0)), _mm256_mul_pd(xi_262, _mm256_set_pd(2.0, 2.0, 2.0, 2.0)));
          const __m256d xi_132 = _mm256_add_pd(xi_118, xi_248);
          const __m256d xi_133 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_132, xi_14), xi_22), xi_251), xi_255);
          const __m256d xi_135 = _mm256_mul_pd(xi_133, _mm256_set_pd(xi_134, xi_134, xi_134, xi_134));
          const __m256d xi_136 = _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_5_1);
          const __m256d xi_141 = _mm256_mul_pd(xi_252, _mm256_set_pd(2.0, 2.0, 2.0, 2.0));
          const __m256d xi_142 = _mm256_mul_pd(xi_247, _mm256_set_pd(2.0, 2.0, 2.0, 2.0));
          const __m256d xi_143 = _mm256_add_pd(_mm256_mul_pd(xi_250, _mm256_set_pd(2.0, 2.0, 2.0, 2.0)), _mm256_mul_pd(xi_257, _mm256_set_pd(-2.0, -2.0, -2.0, -2.0)));
          const __m256d xi_144 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_141, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_14), xi_142), xi_143), xi_19), xi_4);
          const __m256d xi_146 = _mm256_mul_pd(xi_144, _mm256_set_pd(xi_145, xi_145, xi_145, xi_145));
          const __m256d xi_147 = _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_3_0);
          const __m256d xi_152 = _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_0_1);
          const __m256d xi_166 = _mm256_add_pd(xi_122, xi_256);
          const __m256d xi_167 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_12, xi_166), xi_20), xi_244), xi_262);
          const __m256d xi_168 = _mm256_mul_pd(xi_167, _mm256_set_pd(xi_134, xi_134, xi_134, xi_134));
          const __m256d xi_169 = _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_4_1);
          const __m256d xi_171 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_142, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_13), xi_141), xi_143), xi_3);
          const __m256d xi_172 = _mm256_mul_pd(xi_171, _mm256_set_pd(xi_145, xi_145, xi_145, xi_145));
          const __m256d xi_173 = _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_4_0);
          const __m256d xi_178 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_119, xi_23), xi_246), xi_255), xi_265);
          const __m256d xi_179 = _mm256_mul_pd(xi_178, _mm256_set_pd(xi_134, xi_134, xi_134, xi_134));
          const __m256d xi_180 = _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_5_0);
          const __m256d xi_182 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_128, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_129, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_130), xi_24), xi_5);
          const __m256d xi_183 = _mm256_mul_pd(xi_182, _mm256_set_pd(xi_145, xi_145, xi_145, xi_145));
          const __m256d xi_184 = _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_3_1);
          const __m256d xi_212 = _mm256_mul_pd(xi_182, _mm256_set_pd(xi_211, xi_211, xi_211, xi_211));
          const __m256d xi_214 = _mm256_mul_pd(xi_178, _mm256_set_pd(xi_213, xi_213, xi_213, xi_213));
          const __m256d xi_220 = _mm256_mul_pd(xi_144, _mm256_set_pd(xi_211, xi_211, xi_211, xi_211));
          const __m256d xi_221 = _mm256_mul_pd(xi_133, _mm256_set_pd(xi_213, xi_213, xi_213, xi_213));
          const __m256d xi_235 = _mm256_mul_pd(xi_167, _mm256_set_pd(xi_213, xi_213, xi_213, xi_213));
          const __m256d xi_236 = _mm256_mul_pd(xi_171, _mm256_set_pd(xi_211, xi_211, xi_211, xi_211));
          const __m256d xi_31 = _mm256_mul_pd(xi_30, _mm256_set_pd(rr_0, rr_0, rr_0, rr_0));
          const __m256d xi_44 = _mm256_mul_pd(xi_43, _mm256_set_pd(rr_0, rr_0, rr_0, rr_0));
          const __m256d xi_51 = _mm256_mul_pd(xi_50, _mm256_set_pd(rr_0, rr_0, rr_0, rr_0));
          const __m256d xi_54 = _mm256_mul_pd(xi_261, _mm256_set_pd(xi_53, xi_53, xi_53, xi_53));
          const __m256d xi_59 = _mm256_mul_pd(xi_254, _mm256_set_pd(xi_53, xi_53, xi_53, xi_53));
          const __m256d xi_81 = _mm256_mul_pd(xi_259, _mm256_set_pd(xi_53, xi_53, xi_53, xi_53));
          const __m256d vel0Term = _mm256_add_pd(_mm256_add_pd(xi_247, xi_257), xi_3);
          const __m256d vel1Term = _mm256_add_pd(xi_252, xi_4);
          const __m256d vel2Term = _mm256_add_pd(xi_256, xi_5);
          const __m256d rho = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(vel0Term, vel1Term), vel2Term), xi_250), xi_253), xi_255), xi_258), xi_264), xi_6);
          const __m256d xi_105 = _mm256_mul_pd(rho, _mm256_set_pd(kT, kT, kT, kT));
          const __m256d xi_106 = _mm256_sqrt_pd(_mm256_mul_pd(xi_105, _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_set_pd(omega_even, omega_even, omega_even, omega_even)), _mm256_set_pd(1.0, 1.0, 1.0, 1.0)), _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_set_pd(omega_even, omega_even, omega_even, omega_even)), _mm256_set_pd(1.0, 1.0, 1.0, 1.0)))), _mm256_set_pd(1.0, 1.0, 1.0, 1.0))));
          const __m256d xi_107 = _mm256_mul_pd(_mm256_mul_pd(xi_106, _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_6_0)), _mm256_set_pd(3.7416573867739413, 3.7416573867739413, 3.7416573867739413, 3.7416573867739413));
          const __m256d xi_108 = _mm256_mul_pd(_mm256_mul_pd(xi_106, _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_7_0)), _mm256_set_pd(5.4772255750516612, 5.4772255750516612, 5.4772255750516612, 5.4772255750516612));
          const __m256d xi_110 = _mm256_mul_pd(_mm256_mul_pd(_mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_2_1), _mm256_set_pd(xi_109, xi_109, xi_109, xi_109)), _mm256_sqrt_pd(_mm256_mul_pd(xi_105, _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk, omega_bulk)), _mm256_set_pd(1.0, 1.0, 1.0, 1.0)), _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk, omega_bulk)), _mm256_set_pd(1.0, 1.0, 1.0, 1.0)))), _mm256_set_pd(1.0, 1.0, 1.0, 1.0)))));
          const __m256d xi_111 = _mm256_mul_pd(_mm256_mul_pd(xi_106, _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_6_1)), _mm256_set_pd(8.3666002653407556, 8.3666002653407556, 8.3666002653407556, 8.3666002653407556));
          const __m256d xi_137 = _mm256_sqrt_pd(_mm256_mul_pd(xi_105, _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_set_pd(omega_odd, omega_odd, omega_odd, omega_odd)), _mm256_set_pd(1.0, 1.0, 1.0, 1.0)), _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_set_pd(omega_odd, omega_odd, omega_odd, omega_odd)), _mm256_set_pd(1.0, 1.0, 1.0, 1.0)))), _mm256_set_pd(1.0, 1.0, 1.0, 1.0))));
          const __m256d xi_138 = _mm256_mul_pd(xi_137, _mm256_set_pd(1.4142135623730951, 1.4142135623730951, 1.4142135623730951, 1.4142135623730951));
          const __m256d xi_139 = _mm256_mul_pd(xi_138, _mm256_set_pd(0.5, 0.5, 0.5, 0.5));
          const __m256d xi_140 = _mm256_mul_pd(xi_136, xi_139);
          const __m256d xi_148 = _mm256_mul_pd(xi_137, _mm256_set_pd(xi_109, xi_109, xi_109, xi_109));
          const __m256d xi_149 = _mm256_mul_pd(xi_148, _mm256_set_pd(0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666666));
          const __m256d xi_150 = _mm256_mul_pd(xi_147, xi_149);
          const __m256d xi_151 = _mm256_add_pd(_mm256_mul_pd(xi_146, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_150, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)));
          const __m256d xi_153 = _mm256_sqrt_pd(_mm256_mul_pd(xi_105, _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_set_pd(omega_shear, omega_shear, omega_shear, omega_shear)), _mm256_set_pd(1.0, 1.0, 1.0, 1.0)), _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_set_pd(omega_shear, omega_shear, omega_shear, omega_shear)), _mm256_set_pd(1.0, 1.0, 1.0, 1.0)))), _mm256_set_pd(1.0, 1.0, 1.0, 1.0))));
          const __m256d xi_154 = _mm256_mul_pd(xi_153, _mm256_set_pd(0.5, 0.5, 0.5, 0.5));
          const __m256d xi_155 = _mm256_mul_pd(xi_152, xi_154);
          const __m256d xi_161 = _mm256_mul_pd(_mm256_mul_pd(xi_153, _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_0_0)), _mm256_set_pd(1.7320508075688772, 1.7320508075688772, 1.7320508075688772, 1.7320508075688772));
          const __m256d xi_165 = _mm256_add_pd(xi_146, xi_150);
          const __m256d xi_170 = _mm256_mul_pd(xi_139, xi_169);
          const __m256d xi_174 = _mm256_mul_pd(xi_149, xi_173);
          const __m256d xi_175 = _mm256_add_pd(xi_172, xi_174);
          const __m256d xi_177 = _mm256_add_pd(_mm256_mul_pd(xi_172, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_174, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)));
          const __m256d xi_181 = _mm256_mul_pd(xi_139, xi_180);
          const __m256d xi_185 = _mm256_mul_pd(xi_149, xi_184);
          const __m256d xi_186 = _mm256_add_pd(_mm256_mul_pd(xi_183, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_185, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)));
          const __m256d xi_188 = _mm256_add_pd(xi_183, xi_185);
          const __m256d xi_189 = _mm256_mul_pd(_mm256_mul_pd(xi_152, xi_153), _mm256_set_pd(0.25, 0.25, 0.25, 0.25));
          const __m256d xi_192 = _mm256_mul_pd(xi_107, _mm256_set_pd(0.083333333333333329, 0.083333333333333329, 0.083333333333333329, 0.083333333333333329));
          const __m256d xi_196 = _mm256_mul_pd(xi_154, _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_1_0));
          const __m256d xi_203 = _mm256_mul_pd(xi_154, _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_2_0));
          const __m256d xi_207 = _mm256_mul_pd(xi_111, _mm256_set_pd(-0.014285714285714285, -0.014285714285714285, -0.014285714285714285, -0.014285714285714285));
          const __m256d xi_208 = _mm256_mul_pd(xi_108, _mm256_set_pd(0.050000000000000003, 0.050000000000000003, 0.050000000000000003, 0.050000000000000003));
          const __m256d xi_215 = _mm256_mul_pd(xi_148, _mm256_set_pd(0.083333333333333329, 0.083333333333333329, 0.083333333333333329, 0.083333333333333329));
          const __m256d xi_216 = _mm256_mul_pd(xi_184, xi_215);
          const __m256d xi_217 = _mm256_mul_pd(xi_138, _mm256_set_pd(0.25, 0.25, 0.25, 0.25));
          const __m256d xi_218 = _mm256_mul_pd(xi_180, xi_217);
          const __m256d xi_219 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_212, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_216, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_214), xi_218);
          const __m256d xi_222 = _mm256_mul_pd(xi_147, xi_215);
          const __m256d xi_223 = _mm256_mul_pd(xi_136, xi_217);
          const __m256d xi_224 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_220, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_222, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_221), xi_223);
          const __m256d xi_225 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_221, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_223, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_220), xi_222);
          const __m256d xi_227 = _mm256_mul_pd(xi_189, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_230 = _mm256_mul_pd(xi_111, _mm256_set_pd(0.035714285714285712, 0.035714285714285712, 0.035714285714285712, 0.035714285714285712));
          const __m256d xi_232 = _mm256_mul_pd(xi_154, _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_1_1));
          const __m256d xi_237 = _mm256_mul_pd(xi_169, xi_217);
          const __m256d xi_238 = _mm256_mul_pd(xi_173, xi_215);
          const __m256d xi_239 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_235, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_237, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_236), xi_238);
          const __m256d xi_241 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_236, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_238, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_235), xi_237);
          const __m256d xi_242 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_214, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_218, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_212), xi_216);
          const __m256d xi_0 = _mm256_div_pd(_mm256_set_pd(1.0, 1.0, 1.0, 1.0), rho);
          const __m256d xi_7 = _mm256_mul_pd(xi_0, _mm256_set_pd(0.5, 0.5, 0.5, 0.5));
          const __m256d u_0 = _mm256_add_pd(_mm256_mul_pd(xi_0, _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(vel0Term, xi_13), xi_8), xi_9)), _mm256_mul_pd(xi_261, xi_7));
          const __m256d xi_25 = _mm256_mul_pd(u_0, xi_261);
          const __m256d xi_37 = _mm256_mul_pd(xi_25, _mm256_set_pd(0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666666));
          const __m256d xi_38 = _mm256_mul_pd(xi_25, _mm256_set_pd(0.083333333333333329, 0.083333333333333329, 0.083333333333333329, 0.083333333333333329));
          const __m256d xi_39 = _mm256_mul_pd(xi_38, _mm256_set_pd(omega_shear, omega_shear, omega_shear, omega_shear));
          const __m256d xi_40 = _mm256_add_pd(_mm256_mul_pd(xi_37, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_39);
          const __m256d xi_56 = _mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_25, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_set_pd(xi_55, xi_55, xi_55, xi_55)), xi_37);
          const __m256d xi_57 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_43, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_54), xi_56);
          const __m256d xi_61 = _mm256_mul_pd(_mm256_mul_pd(xi_25, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_set_pd(xi_60, xi_60, xi_60, xi_60));
          const __m256d xi_68 = _mm256_mul_pd(u_0, xi_67);
          const __m256d xi_73 = _mm256_mul_pd(u_0, xi_72);
          const __m256d xi_77 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_54, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_43), xi_56);
          const __m256d xi_84 = _mm256_mul_pd(xi_38, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_95 = _mm256_mul_pd(u_0, xi_259);
          const __m256d xi_96 = _mm256_mul_pd(xi_95, _mm256_set_pd(0.25, 0.25, 0.25, 0.25));
          const __m256d xi_99 = _mm256_mul_pd(xi_95, _mm256_set_pd(xi_71, xi_71, xi_71, xi_71));
          const __m256d xi_113 = _mm256_mul_pd(rho, _mm256_mul_pd(u_0, u_0));
          const __m256d u_1 = _mm256_add_pd(_mm256_mul_pd(xi_0, _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(vel1Term, xi_16), xi_19), xi_257), xi_8)), _mm256_mul_pd(xi_254, xi_7));
          const __m256d xi_26 = _mm256_mul_pd(u_1, xi_254);
          const __m256d xi_32 = _mm256_mul_pd(xi_26, _mm256_set_pd(0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666666));
          const __m256d xi_45 = _mm256_mul_pd(xi_26, _mm256_set_pd(0.083333333333333329, 0.083333333333333329, 0.083333333333333329, 0.083333333333333329));
          const __m256d xi_46 = _mm256_mul_pd(xi_45, _mm256_set_pd(omega_shear, omega_shear, omega_shear, omega_shear));
          const __m256d xi_47 = _mm256_add_pd(_mm256_mul_pd(xi_32, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_46);
          const __m256d xi_62 = _mm256_mul_pd(_mm256_mul_pd(xi_26, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_set_pd(xi_60, xi_60, xi_60, xi_60));
          const __m256d xi_69 = _mm256_mul_pd(u_1, _mm256_set_pd(0.25, 0.25, 0.25, 0.25));
          const __m256d xi_70 = _mm256_mul_pd(xi_261, xi_69);
          const __m256d xi_74 = _mm256_mul_pd(u_1, _mm256_set_pd(xi_71, xi_71, xi_71, xi_71));
          const __m256d xi_75 = _mm256_mul_pd(xi_261, xi_74);
          const __m256d xi_76 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_68, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_70, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_73), xi_75);
          const __m256d xi_78 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_73, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_75, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_68), xi_70);
          const __m256d xi_86 = _mm256_mul_pd(xi_259, xi_69);
          const __m256d xi_88 = _mm256_mul_pd(xi_259, xi_74);
          const __m256d xi_93 = _mm256_mul_pd(xi_45, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_112 = _mm256_mul_pd(rho, _mm256_mul_pd(u_1, u_1));
          const __m256d xi_121 = _mm256_add_pd(_mm256_add_pd(xi_112, xi_120), xi_9);
          const __m256d xi_197 = _mm256_mul_pd(rho, u_1);
          const __m256d xi_199 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(u_0, xi_197), xi_120), xi_247), xi_252), _mm256_set_pd(xi_198, xi_198, xi_198, xi_198));
          const __m256d xi_200 = _mm256_add_pd(_mm256_mul_pd(xi_196, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_199, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)));
          const __m256d xi_201 = _mm256_add_pd(xi_196, xi_199);
          const __m256d u_2 = _mm256_add_pd(_mm256_mul_pd(xi_0, _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(vel2Term, xi_21), xi_24), xi_260)), _mm256_mul_pd(xi_259, xi_7));
          const __m256d xi_27 = _mm256_mul_pd(u_2, xi_259);
          const __m256d xi_33 = _mm256_mul_pd(xi_27, _mm256_set_pd(0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666666));
          const __m256d xi_34 = _mm256_mul_pd(xi_27, _mm256_set_pd(0.083333333333333329, 0.083333333333333329, 0.083333333333333329, 0.083333333333333329));
          const __m256d xi_35 = _mm256_mul_pd(xi_34, _mm256_set_pd(omega_shear, omega_shear, omega_shear, omega_shear));
          const __m256d xi_36 = _mm256_add_pd(_mm256_mul_pd(xi_33, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_35);
          const __m256d xi_41 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_26, _mm256_set_pd(0.33333333333333331, 0.33333333333333331, 0.33333333333333331, 0.33333333333333331)), _mm256_mul_pd(_mm256_mul_pd(xi_32, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_set_pd(omega_shear, omega_shear, omega_shear, omega_shear))), xi_36), xi_40);
          const __m256d xi_48 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_25, _mm256_set_pd(0.33333333333333331, 0.33333333333333331, 0.33333333333333331, 0.33333333333333331)), _mm256_mul_pd(_mm256_mul_pd(xi_37, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_set_pd(omega_shear, omega_shear, omega_shear, omega_shear))), xi_36), xi_47);
          const __m256d xi_52 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_27, _mm256_set_pd(0.33333333333333331, 0.33333333333333331, 0.33333333333333331, 0.33333333333333331)), _mm256_mul_pd(_mm256_mul_pd(xi_33, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_set_pd(omega_shear, omega_shear, omega_shear, omega_shear))), xi_40), xi_47);
          const __m256d xi_58 = _mm256_mul_pd(xi_34, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_63 = _mm256_mul_pd(_mm256_mul_pd(xi_27, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_set_pd(xi_60, xi_60, xi_60, xi_60));
          const __m256d xi_64 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_26, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_set_pd(xi_55, xi_55, xi_55, xi_55)), xi_32), xi_61), xi_62), xi_63);
          const __m256d xi_65 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_59, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_30), xi_64);
          const __m256d xi_66 = _mm256_add_pd(_mm256_add_pd(xi_35, xi_58), xi_65);
          const __m256d xi_79 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_30, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_59), xi_64);
          const __m256d xi_80 = _mm256_add_pd(_mm256_add_pd(xi_35, xi_58), xi_79);
          const __m256d xi_82 = _mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_27, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_set_pd(xi_55, xi_55, xi_55, xi_55)), xi_33);
          const __m256d xi_83 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_81, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_50), xi_82);
          const __m256d xi_85 = _mm256_add_pd(_mm256_add_pd(xi_39, xi_65), xi_84);
          const __m256d xi_87 = _mm256_mul_pd(u_2, xi_67);
          const __m256d xi_89 = _mm256_mul_pd(u_2, xi_72);
          const __m256d xi_90 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_88, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_89, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_86), xi_87);
          const __m256d xi_91 = _mm256_add_pd(_mm256_add_pd(xi_39, xi_79), xi_84);
          const __m256d xi_92 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_86, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_87, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_88), xi_89);
          const __m256d xi_94 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_46, xi_61), xi_62), xi_63), xi_83), xi_93);
          const __m256d xi_97 = _mm256_mul_pd(u_2, xi_261);
          const __m256d xi_98 = _mm256_mul_pd(xi_97, _mm256_set_pd(0.25, 0.25, 0.25, 0.25));
          const __m256d xi_100 = _mm256_mul_pd(xi_97, _mm256_set_pd(xi_71, xi_71, xi_71, xi_71));
          const __m256d xi_101 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_96, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_98, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_100), xi_99);
          const __m256d xi_102 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_100, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_99, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_96), xi_98);
          const __m256d xi_103 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_50, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_81), xi_82);
          const __m256d xi_104 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_103, xi_46), xi_61), xi_62), xi_63), xi_93);
          const __m256d xi_115 = _mm256_mul_pd(rho, _mm256_mul_pd(u_2, u_2));
          const __m256d xi_116 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_245, _mm256_set_pd(3.0, 3.0, 3.0, 3.0)), _mm256_mul_pd(xi_265, _mm256_set_pd(3.0, 3.0, 3.0, 3.0))), _mm256_mul_pd(xi_115, _mm256_set_pd(0.66666666666666663, 0.66666666666666663, 0.66666666666666663, 0.66666666666666663))), xi_114);
          const __m256d xi_117 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_251, _mm256_set_pd(3.0, 3.0, 3.0, 3.0)), _mm256_mul_pd(xi_264, _mm256_set_pd(3.0, 3.0, 3.0, 3.0))), _mm256_mul_pd(xi_112, _mm256_set_pd(0.66666666666666663, 0.66666666666666663, 0.66666666666666663, 0.66666666666666663))), _mm256_mul_pd(xi_113, _mm256_set_pd(1.6666666666666667, 1.6666666666666667, 1.6666666666666667, 1.6666666666666667))), _mm256_mul_pd(xi_246, _mm256_set_pd(-3.0, -3.0, -3.0, -3.0))), _mm256_mul_pd(xi_248, _mm256_set_pd(-3.0, -3.0, -3.0, -3.0))), _mm256_mul_pd(xi_255, _mm256_set_pd(-3.0, -3.0, -3.0, -3.0))), _mm256_mul_pd(xi_263, _mm256_set_pd(-3.0, -3.0, -3.0, -3.0))), xi_116), _mm256_set_pd(omega_even, omega_even, omega_even, omega_even));
          const __m256d xi_124 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_113, xi_115), xi_119), xi_121), xi_123), xi_17), xi_22), xi_253), _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk, omega_bulk));
          const __m256d xi_127 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_112, _mm256_set_pd(2.3333333333333335, 2.3333333333333335, 2.3333333333333335, 2.3333333333333335)), _mm256_mul_pd(xi_251, _mm256_set_pd(-2.0, -2.0, -2.0, -2.0))), _mm256_mul_pd(xi_264, _mm256_set_pd(-2.0, -2.0, -2.0, -2.0))), _mm256_mul_pd(xi_249, _mm256_set_pd(-5.0, -5.0, -5.0, -5.0))), _mm256_mul_pd(xi_256, _mm256_set_pd(-5.0, -5.0, -5.0, -5.0))), _mm256_mul_pd(xi_260, _mm256_set_pd(-5.0, -5.0, -5.0, -5.0))), _mm256_mul_pd(xi_262, _mm256_set_pd(-5.0, -5.0, -5.0, -5.0))), xi_116), xi_126), _mm256_set_pd(omega_even, omega_even, omega_even, omega_even));
          const __m256d xi_131 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_115, _mm256_set_pd(3.0, 3.0, 3.0, 3.0)), _mm256_mul_pd(xi_251, _mm256_set_pd(5.0, 5.0, 5.0, 5.0))), _mm256_mul_pd(xi_264, _mm256_set_pd(5.0, 5.0, 5.0, 5.0))), _mm256_mul_pd(xi_245, _mm256_set_pd(-4.0, -4.0, -4.0, -4.0))), _mm256_mul_pd(xi_265, _mm256_set_pd(-4.0, -4.0, -4.0, -4.0))), _mm256_mul_pd(xi_247, _mm256_set_pd(-7.0, -7.0, -7.0, -7.0))), _mm256_mul_pd(xi_250, _mm256_set_pd(-7.0, -7.0, -7.0, -7.0))), _mm256_mul_pd(xi_252, _mm256_set_pd(-7.0, -7.0, -7.0, -7.0))), _mm256_mul_pd(xi_257, _mm256_set_pd(-7.0, -7.0, -7.0, -7.0))), xi_114), xi_126), xi_128), xi_129), xi_130), _mm256_set_pd(omega_even, omega_even, omega_even, omega_even));
          const __m256d xi_156 = _mm256_add_pd(_mm256_mul_pd(xi_115, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_265);
          const __m256d xi_157 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_251, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_121), xi_156), xi_16), xi_2), xi_256), xi_6), _mm256_set_pd(omega_shear, omega_shear, omega_shear, omega_shear));
          const __m256d xi_158 = _mm256_mul_pd(xi_157, _mm256_set_pd(0.125, 0.125, 0.125, 0.125));
          const __m256d xi_159 = _mm256_add_pd(_mm256_mul_pd(xi_131, _mm256_set_pd(-0.01984126984126984, -0.01984126984126984, -0.01984126984126984, -0.01984126984126984)), _mm256_mul_pd(xi_107, _mm256_set_pd(-0.11904761904761904, -0.11904761904761904, -0.11904761904761904, -0.11904761904761904)));
          const __m256d xi_160 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_112, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_113, _mm256_set_pd(2.0, 2.0, 2.0, 2.0))), _mm256_mul_pd(xi_244, _mm256_set_pd(-2.0, -2.0, -2.0, -2.0))), _mm256_mul_pd(xi_258, _mm256_set_pd(-2.0, -2.0, -2.0, -2.0))), xi_120), xi_123), xi_125), xi_156), xi_245), xi_251), xi_264), xi_9), _mm256_set_pd(omega_shear, omega_shear, omega_shear, omega_shear));
          const __m256d xi_162 = _mm256_add_pd(_mm256_mul_pd(xi_160, _mm256_set_pd(-0.041666666666666664, -0.041666666666666664, -0.041666666666666664, -0.041666666666666664)), _mm256_mul_pd(xi_161, _mm256_set_pd(-0.16666666666666666, -0.16666666666666666, -0.16666666666666666, -0.16666666666666666)));
          const __m256d xi_163 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_117, _mm256_set_pd(-0.050000000000000003, -0.050000000000000003, -0.050000000000000003, -0.050000000000000003)), _mm256_mul_pd(xi_108, _mm256_set_pd(-0.10000000000000001, -0.10000000000000001, -0.10000000000000001, -0.10000000000000001))), xi_162);
          const __m256d xi_164 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_127, _mm256_set_pd(0.014285714285714285, 0.014285714285714285, 0.014285714285714285, 0.014285714285714285)), _mm256_mul_pd(xi_111, _mm256_set_pd(0.028571428571428571, 0.028571428571428571, 0.028571428571428571, 0.028571428571428571))), xi_155), xi_158), xi_159), xi_163);
          const __m256d xi_176 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_160, _mm256_set_pd(0.083333333333333329, 0.083333333333333329, 0.083333333333333329, 0.083333333333333329)), _mm256_mul_pd(xi_161, _mm256_set_pd(0.33333333333333331, 0.33333333333333331, 0.33333333333333331, 0.33333333333333331))), _mm256_mul_pd(xi_127, _mm256_set_pd(-0.035714285714285712, -0.035714285714285712, -0.035714285714285712, -0.035714285714285712))), _mm256_mul_pd(xi_111, _mm256_set_pd(-0.071428571428571425, -0.071428571428571425, -0.071428571428571425, -0.071428571428571425))), xi_159);
          const __m256d xi_187 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_155, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_158, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_131, _mm256_set_pd(0.015873015873015872, 0.015873015873015872, 0.015873015873015872, 0.015873015873015872))), _mm256_mul_pd(xi_107, _mm256_set_pd(0.095238095238095233, 0.095238095238095233, 0.095238095238095233, 0.095238095238095233))), _mm256_mul_pd(xi_127, _mm256_set_pd(-0.021428571428571429, -0.021428571428571429, -0.021428571428571429, -0.021428571428571429))), _mm256_mul_pd(xi_111, _mm256_set_pd(-0.042857142857142858, -0.042857142857142858, -0.042857142857142858, -0.042857142857142858))), xi_163);
          const __m256d xi_190 = _mm256_mul_pd(xi_157, _mm256_set_pd(0.0625, 0.0625, 0.0625, 0.0625));
          const __m256d xi_191 = _mm256_mul_pd(xi_131, _mm256_set_pd(0.013888888888888888, 0.013888888888888888, 0.013888888888888888, 0.013888888888888888));
          const __m256d xi_193 = _mm256_add_pd(_mm256_mul_pd(xi_124, _mm256_set_pd(0.041666666666666664, 0.041666666666666664, 0.041666666666666664, 0.041666666666666664)), _mm256_mul_pd(xi_110, _mm256_set_pd(0.083333333333333329, 0.083333333333333329, 0.083333333333333329, 0.083333333333333329)));
          const __m256d xi_194 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_160, _mm256_set_pd(0.020833333333333332, 0.020833333333333332, 0.020833333333333332, 0.020833333333333332)), _mm256_mul_pd(xi_161, _mm256_set_pd(0.083333333333333329, 0.083333333333333329, 0.083333333333333329, 0.083333333333333329))), xi_193);
          const __m256d xi_195 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_165, xi_189), xi_190), xi_191), xi_192), xi_194);
          const __m256d xi_202 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_151, xi_189), xi_190), xi_191), xi_192), xi_194);
          const __m256d xi_204 = _mm256_mul_pd(xi_127, _mm256_set_pd(-0.0071428571428571426, -0.0071428571428571426, -0.0071428571428571426, -0.0071428571428571426));
          const __m256d xi_205 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(u_2, xi_197), xi_132), xi_17), xi_246), _mm256_set_pd(xi_198, xi_198, xi_198, xi_198));
          const __m256d xi_206 = _mm256_mul_pd(xi_117, _mm256_set_pd(0.025000000000000001, 0.025000000000000001, 0.025000000000000001, 0.025000000000000001));
          const __m256d xi_209 = _mm256_add_pd(_mm256_mul_pd(xi_131, _mm256_set_pd(-0.003968253968253968, -0.003968253968253968, -0.003968253968253968, -0.003968253968253968)), _mm256_mul_pd(xi_107, _mm256_set_pd(-0.023809523809523808, -0.023809523809523808, -0.023809523809523808, -0.023809523809523808)));
          const __m256d xi_210 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_162, xi_193), xi_203), xi_204), xi_205), xi_206), xi_207), xi_208), xi_209);
          const __m256d xi_226 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_203, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_205, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_162), xi_193), xi_204), xi_206), xi_207), xi_208), xi_209);
          const __m256d xi_228 = _mm256_mul_pd(xi_190, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
          const __m256d xi_229 = _mm256_mul_pd(xi_127, _mm256_set_pd(0.017857142857142856, 0.017857142857142856, 0.017857142857142856, 0.017857142857142856));
          const __m256d xi_231 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_188, xi_194), xi_209), xi_227), xi_228), xi_229), xi_230);
          const __m256d xi_233 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(rho, u_0), u_2), xi_10), xi_166), xi_249), _mm256_set_pd(xi_198, xi_198, xi_198, xi_198));
          const __m256d xi_234 = _mm256_add_pd(_mm256_mul_pd(xi_232, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_233, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)));
          const __m256d xi_240 = _mm256_add_pd(xi_232, xi_233);
          const __m256d xi_243 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_186, xi_194), xi_209), xi_227), xi_228), xi_229), xi_230);
          const __m256d forceTerm_0 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_25, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_26, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_27, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_25, _mm256_set_pd(xi_28, xi_28, xi_28, xi_28))), _mm256_mul_pd(xi_26, _mm256_set_pd(xi_28, xi_28, xi_28, xi_28))), _mm256_mul_pd(xi_27, _mm256_set_pd(xi_28, xi_28, xi_28, xi_28)));
          const __m256d forceTerm_1 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_31, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_29), xi_41);
          const __m256d forceTerm_2 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_29, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_31), xi_41);
          const __m256d forceTerm_3 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_42, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_44), xi_48);
          const __m256d forceTerm_4 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_44, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_42), xi_48);
          const __m256d forceTerm_5 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_51, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_49), xi_52);
          const __m256d forceTerm_6 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_49, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_51), xi_52);
          const __m256d forceTerm_7 = _mm256_add_pd(_mm256_add_pd(xi_57, xi_66), xi_76);
          const __m256d forceTerm_8 = _mm256_add_pd(_mm256_add_pd(xi_66, xi_77), xi_78);
          const __m256d forceTerm_9 = _mm256_add_pd(_mm256_add_pd(xi_57, xi_78), xi_80);
          const __m256d forceTerm_10 = _mm256_add_pd(_mm256_add_pd(xi_76, xi_77), xi_80);
          const __m256d forceTerm_11 = _mm256_add_pd(_mm256_add_pd(xi_83, xi_85), xi_90);
          const __m256d forceTerm_12 = _mm256_add_pd(_mm256_add_pd(xi_83, xi_91), xi_92);
          const __m256d forceTerm_13 = _mm256_add_pd(_mm256_add_pd(xi_101, xi_57), xi_94);
          const __m256d forceTerm_14 = _mm256_add_pd(_mm256_add_pd(xi_102, xi_77), xi_94);
          const __m256d forceTerm_15 = _mm256_add_pd(_mm256_add_pd(xi_103, xi_85), xi_92);
          const __m256d forceTerm_16 = _mm256_add_pd(_mm256_add_pd(xi_103, xi_90), xi_91);
          const __m256d forceTerm_17 = _mm256_add_pd(_mm256_add_pd(xi_102, xi_104), xi_57);
          const __m256d forceTerm_18 = _mm256_add_pd(_mm256_add_pd(xi_101, xi_104), xi_77);
          _mm256_store_pd(&_data_pdfs_20_30_10[ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_110, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_131, _mm256_set_pd(0.023809523809523808, 0.023809523809523808, 0.023809523809523808, 0.023809523809523808))), _mm256_mul_pd(xi_107, _mm256_set_pd(0.14285714285714285, 0.14285714285714285, 0.14285714285714285, 0.14285714285714285))), _mm256_mul_pd(xi_127, _mm256_set_pd(0.042857142857142858, 0.042857142857142858, 0.042857142857142858, 0.042857142857142858))), _mm256_mul_pd(xi_111, _mm256_set_pd(0.085714285714285715, 0.085714285714285715, 0.085714285714285715, 0.085714285714285715))), _mm256_mul_pd(xi_117, _mm256_set_pd(0.10000000000000001, 0.10000000000000001, 0.10000000000000001, 0.10000000000000001))), _mm256_mul_pd(xi_108, _mm256_set_pd(0.20000000000000001, 0.20000000000000001, 0.20000000000000001, 0.20000000000000001))), _mm256_mul_pd(xi_124, _mm256_set_pd(-0.5, -0.5, -0.5, -0.5))), forceTerm_0), xi_253));
          _mm256_store_pd(&_data_pdfs_20_31_10[ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_135, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_140, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), forceTerm_1), xi_151), xi_164), xi_251));
          _mm256_store_pd(&_data_pdfs_20_32_10[ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_2, xi_135), xi_140), xi_164), xi_165), xi_264));
          _mm256_store_pd(&_data_pdfs_20_33_10[ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_3, xi_168), xi_170), xi_175), xi_176), xi_258));
          _mm256_store_pd(&_data_pdfs_20_34_10[ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_168, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_170, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), forceTerm_4), xi_176), xi_177), xi_244));
          _mm256_store_pd(&_data_pdfs_20_35_10[ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_179, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_181, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), forceTerm_5), xi_186), xi_187), xi_265));
          _mm256_store_pd(&_data_pdfs_20_36_10[ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_6, xi_179), xi_181), xi_187), xi_188), xi_245));
          _mm256_store_pd(&_data_pdfs_20_37_10[ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_7, xi_177), xi_195), xi_200), xi_252));
          _mm256_store_pd(&_data_pdfs_20_38_10[ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_8, xi_175), xi_195), xi_201), xi_257));
          _mm256_store_pd(&_data_pdfs_20_39_10[ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_9, xi_177), xi_201), xi_202), xi_250));
          _mm256_store_pd(&_data_pdfs_20_310_10[ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_10, xi_175), xi_200), xi_202), xi_247));
          _mm256_store_pd(&_data_pdfs_20_311_10[ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_11, xi_210), xi_219), xi_224), xi_263));
          _mm256_store_pd(&_data_pdfs_20_312_10[ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_12, xi_219), xi_225), xi_226), xi_248));
          _mm256_store_pd(&_data_pdfs_20_313_10[ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_13, xi_231), xi_234), xi_239), xi_256));
          _mm256_store_pd(&_data_pdfs_20_314_10[ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_14, xi_231), xi_240), xi_241), xi_260));
          _mm256_store_pd(&_data_pdfs_20_315_10[ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_15, xi_224), xi_226), xi_242), xi_246));
          _mm256_store_pd(&_data_pdfs_20_316_10[ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_16, xi_210), xi_225), xi_242), xi_255));
          _mm256_store_pd(&_data_pdfs_20_317_10[ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_17, xi_239), xi_240), xi_243), xi_262));
          _mm256_store_pd(&_data_pdfs_20_318_10[ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_18, xi_234), xi_241), xi_243), xi_249));
        }
        for (int64_t ctr_0 = (int64_t)((_size_force_0) / (4)) * (4); ctr_0 < _size_force_0; ctr_0 += 1) {
          const double xi_244 = _data_pdfs_20_34_10[ctr_0];
          const double xi_245 = _data_pdfs_20_36_10[ctr_0];
          const double xi_246 = _data_pdfs_20_315_10[ctr_0];
          const double xi_247 = _data_pdfs_20_310_10[ctr_0];
          const double xi_248 = _data_pdfs_20_312_10[ctr_0];
          const double xi_249 = _data_pdfs_20_318_10[ctr_0];
          const double xi_250 = _data_pdfs_20_39_10[ctr_0];
          const double xi_251 = _data_pdfs_20_31_10[ctr_0];
          const double xi_252 = _data_pdfs_20_37_10[ctr_0];
          const double xi_253 = _data_pdfs_20_30_10[ctr_0];
          const double xi_254 = _data_force_20_31_10[ctr_0];
          const double xi_255 = _data_pdfs_20_316_10[ctr_0];
          const double xi_256 = _data_pdfs_20_313_10[ctr_0];
          const double xi_257 = _data_pdfs_20_38_10[ctr_0];
          const double xi_258 = _data_pdfs_20_33_10[ctr_0];
          const double xi_259 = _data_force_20_32_10[ctr_0];
          const double xi_260 = _data_pdfs_20_314_10[ctr_0];
          const double xi_261 = _data_force_20_30_10[ctr_0];
          const double xi_262 = _data_pdfs_20_317_10[ctr_0];
          const double xi_263 = _data_pdfs_20_311_10[ctr_0];
          const double xi_264 = _data_pdfs_20_32_10[ctr_0];
          const double xi_265 = _data_pdfs_20_35_10[ctr_0];

          double random_7_0{};
          double random_7_1{};
          if (kT > 0.) {
            philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 7, seed, random_7_0, random_7_1);
          }

          double random_6_0{};
          double random_6_1{};
          if (kT > 0.) {
            philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 6, seed, random_6_0, random_6_1);
          }

          double random_5_0{};
          double random_5_1{};
          if (kT > 0.) {
            philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 5, seed, random_5_0, random_5_1);
          }

          double random_4_0{};
          double random_4_1{};
          if (kT > 0.) {
            philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 4, seed, random_4_0, random_4_1);
          }

          double random_3_0{};
          double random_3_1{};
          if (kT > 0.) {
            philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 3, seed, random_3_0, random_3_1);
          }

          double random_2_0{};
          double random_2_1{};
          if (kT > 0.) {
            philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 2, seed, random_2_0, random_2_1);
          }

          double random_1_0{};
          double random_1_1{};
          if (kT > 0.) {
            philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 1, seed, random_1_0, random_1_1);
          }

          double random_0_0{};
          double random_0_1{};
          if (kT > 0.) {
            philox_double2(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1, block_offset_2 + ctr_2, 0, seed, random_0_0, random_0_1);
          }
          const double xi_2 = xi_249 + xi_260;
          const double xi_3 = xi_2 + xi_244;
          const double xi_4 = xi_246 + xi_251 + xi_263;
          const double xi_5 = xi_248 + xi_265;
          const double xi_6 = xi_245 + xi_262;
          const double xi_8 = xi_250 * -1.0;
          const double xi_9 = xi_252 * -1.0;
          const double xi_10 = xi_262 * -1.0;
          const double xi_11 = xi_256 * -1.0;
          const double xi_12 = xi_258 * -1.0;
          const double xi_13 = xi_10 + xi_11 + xi_12;
          const double xi_14 = xi_264 * -1.0;
          const double xi_15 = xi_247 * -1.0;
          const double xi_16 = xi_14 + xi_15;
          const double xi_17 = xi_255 * -1.0;
          const double xi_18 = xi_248 * -1.0;
          const double xi_19 = xi_17 + xi_18;
          const double xi_20 = xi_249 * -1.0;
          const double xi_21 = xi_10 + xi_20;
          const double xi_22 = xi_246 * -1.0;
          const double xi_23 = xi_245 * -1.0;
          const double xi_24 = xi_17 + xi_22 + xi_23 + xi_263;
          const double xi_29 = xi_254 * 0.16666666666666666;
          const double xi_30 = xi_254 * 0.083333333333333329;
          const double xi_42 = xi_261 * 0.16666666666666666;
          const double xi_43 = xi_261 * 0.083333333333333329;
          const double xi_49 = xi_259 * 0.16666666666666666;
          const double xi_50 = xi_259 * 0.083333333333333329;
          const double xi_67 = xi_254 * 0.25;
          const double xi_72 = xi_254 * xi_71;
          const double xi_114 = xi_253 * -1.0;
          const double xi_118 = xi_263 * -1.0;
          const double xi_119 = xi_118 + xi_18;
          const double xi_120 = xi_257 * -1.0 + xi_8;
          const double xi_122 = xi_260 * -1.0;
          const double xi_123 = xi_11 + xi_122 + xi_15 + xi_21;
          const double xi_125 = xi_246 * 2.0 + xi_248 * 2.0 + xi_255 * 2.0 + xi_263 * 2.0;
          const double xi_126 = xi_125 + xi_244 * 5.0 + xi_258 * 5.0;
          const double xi_128 = xi_256 * 2.0;
          const double xi_129 = xi_260 * 2.0;
          const double xi_130 = xi_249 * 2.0 + xi_262 * 2.0;
          const double xi_132 = xi_118 + xi_248;
          const double xi_133 = xi_132 + xi_14 + xi_22 + xi_251 + xi_255;
          const double xi_135 = xi_133 * xi_134;
          const double xi_136 = random_5_1 - 0.5;
          const double xi_141 = xi_252 * 2.0;
          const double xi_142 = xi_247 * 2.0;
          const double xi_143 = xi_250 * 2.0 + xi_257 * -2.0;
          const double xi_144 = xi_14 + xi_141 * -1.0 + xi_142 + xi_143 + xi_19 + xi_4;
          const double xi_146 = xi_144 * xi_145;
          const double xi_147 = random_3_0 - 0.5;
          const double xi_152 = random_0_1 - 0.5;
          const double xi_166 = xi_122 + xi_256;
          const double xi_167 = xi_12 + xi_166 + xi_20 + xi_244 + xi_262;
          const double xi_168 = xi_134 * xi_167;
          const double xi_169 = random_4_1 - 0.5;
          const double xi_171 = xi_13 + xi_141 + xi_142 * -1.0 + xi_143 + xi_3;
          const double xi_172 = xi_145 * xi_171;
          const double xi_173 = random_4_0 - 0.5;
          const double xi_178 = xi_119 + xi_23 + xi_246 + xi_255 + xi_265;
          const double xi_179 = xi_134 * xi_178;
          const double xi_180 = random_5_0 - 0.5;
          const double xi_182 = xi_128 * -1.0 + xi_129 * -1.0 + xi_130 + xi_24 + xi_5;
          const double xi_183 = xi_145 * xi_182;
          const double xi_184 = random_3_1 - 0.5;
          const double xi_212 = xi_182 * xi_211;
          const double xi_214 = xi_178 * xi_213;
          const double xi_220 = xi_144 * xi_211;
          const double xi_221 = xi_133 * xi_213;
          const double xi_235 = xi_167 * xi_213;
          const double xi_236 = xi_171 * xi_211;
          const double xi_31 = rr_0 * xi_30;
          const double xi_44 = rr_0 * xi_43;
          const double xi_51 = rr_0 * xi_50;
          const double xi_54 = xi_261 * xi_53;
          const double xi_59 = xi_254 * xi_53;
          const double xi_81 = xi_259 * xi_53;
          const double vel0Term = xi_247 + xi_257 + xi_3;
          const double vel1Term = xi_252 + xi_4;
          const double vel2Term = xi_256 + xi_5;
          const double rho = vel0Term + vel1Term + vel2Term + xi_250 + xi_253 + xi_255 + xi_258 + xi_264 + xi_6;
          const double xi_105 = kT * rho;
          const double xi_106 = pow(xi_105 * (-1.0 * (omega_even * -1.0 + 1.0) * (omega_even * -1.0 + 1.0) + 1.0), 0.5);
          const double xi_107 = xi_106 * (random_6_0 - 0.5) * 3.7416573867739413;
          const double xi_108 = xi_106 * (random_7_0 - 0.5) * 5.4772255750516612;
          const double xi_110 = xi_109 * (random_2_1 - 0.5) * pow(xi_105 * (-1.0 * (omega_bulk * -1.0 + 1.0) * (omega_bulk * -1.0 + 1.0) + 1.0), 0.5);
          const double xi_111 = xi_106 * (random_6_1 - 0.5) * 8.3666002653407556;
          const double xi_137 = pow(xi_105 * (-1.0 * (omega_odd * -1.0 + 1.0) * (omega_odd * -1.0 + 1.0) + 1.0), 0.5);
          const double xi_138 = xi_137 * 1.4142135623730951;
          const double xi_139 = xi_138 * 0.5;
          const double xi_140 = xi_136 * xi_139;
          const double xi_148 = xi_109 * xi_137;
          const double xi_149 = xi_148 * 0.16666666666666666;
          const double xi_150 = xi_147 * xi_149;
          const double xi_151 = xi_146 * -1.0 + xi_150 * -1.0;
          const double xi_153 = pow(xi_105 * (-1.0 * (omega_shear * -1.0 + 1.0) * (omega_shear * -1.0 + 1.0) + 1.0), 0.5);
          const double xi_154 = xi_153 * 0.5;
          const double xi_155 = xi_152 * xi_154;
          const double xi_161 = xi_153 * (random_0_0 - 0.5) * 1.7320508075688772;
          const double xi_165 = xi_146 + xi_150;
          const double xi_170 = xi_139 * xi_169;
          const double xi_174 = xi_149 * xi_173;
          const double xi_175 = xi_172 + xi_174;
          const double xi_177 = xi_172 * -1.0 + xi_174 * -1.0;
          const double xi_181 = xi_139 * xi_180;
          const double xi_185 = xi_149 * xi_184;
          const double xi_186 = xi_183 * -1.0 + xi_185 * -1.0;
          const double xi_188 = xi_183 + xi_185;
          const double xi_189 = xi_152 * xi_153 * 0.25;
          const double xi_192 = xi_107 * 0.083333333333333329;
          const double xi_196 = xi_154 * (random_1_0 - 0.5);
          const double xi_203 = xi_154 * (random_2_0 - 0.5);
          const double xi_207 = xi_111 * -0.014285714285714285;
          const double xi_208 = xi_108 * 0.050000000000000003;
          const double xi_215 = xi_148 * 0.083333333333333329;
          const double xi_216 = xi_184 * xi_215;
          const double xi_217 = xi_138 * 0.25;
          const double xi_218 = xi_180 * xi_217;
          const double xi_219 = xi_212 * -1.0 + xi_214 + xi_216 * -1.0 + xi_218;
          const double xi_222 = xi_147 * xi_215;
          const double xi_223 = xi_136 * xi_217;
          const double xi_224 = xi_220 * -1.0 + xi_221 + xi_222 * -1.0 + xi_223;
          const double xi_225 = xi_220 + xi_221 * -1.0 + xi_222 + xi_223 * -1.0;
          const double xi_227 = xi_189 * -1.0;
          const double xi_230 = xi_111 * 0.035714285714285712;
          const double xi_232 = xi_154 * (random_1_1 - 0.5);
          const double xi_237 = xi_169 * xi_217;
          const double xi_238 = xi_173 * xi_215;
          const double xi_239 = xi_235 * -1.0 + xi_236 + xi_237 * -1.0 + xi_238;
          const double xi_241 = xi_235 + xi_236 * -1.0 + xi_237 + xi_238 * -1.0;
          const double xi_242 = xi_212 + xi_214 * -1.0 + xi_216 + xi_218 * -1.0;
          const double xi_0 = ((1.0) / (rho));
          const double xi_7 = xi_0 * 0.5;
          const double u_0 = xi_0 * (vel0Term + xi_13 + xi_8 + xi_9) + xi_261 * xi_7;
          const double xi_25 = u_0 * xi_261;
          const double xi_37 = xi_25 * 0.16666666666666666;
          const double xi_38 = xi_25 * 0.083333333333333329;
          const double xi_39 = omega_shear * xi_38;
          const double xi_40 = xi_37 * -1.0 + xi_39;
          const double xi_56 = xi_25 * xi_55 * -1.0 + xi_37;
          const double xi_57 = xi_43 * -1.0 + xi_54 + xi_56;
          const double xi_61 = xi_25 * xi_60 * -1.0;
          const double xi_68 = u_0 * xi_67;
          const double xi_73 = u_0 * xi_72;
          const double xi_77 = xi_43 + xi_54 * -1.0 + xi_56;
          const double xi_84 = xi_38 * -1.0;
          const double xi_95 = u_0 * xi_259;
          const double xi_96 = xi_95 * 0.25;
          const double xi_99 = xi_71 * xi_95;
          const double xi_113 = rho * u_0 * u_0;
          const double u_1 = xi_0 * (vel1Term + xi_16 + xi_19 + xi_257 + xi_8) + xi_254 * xi_7;
          const double xi_26 = u_1 * xi_254;
          const double xi_32 = xi_26 * 0.16666666666666666;
          const double xi_45 = xi_26 * 0.083333333333333329;
          const double xi_46 = omega_shear * xi_45;
          const double xi_47 = xi_32 * -1.0 + xi_46;
          const double xi_62 = xi_26 * xi_60 * -1.0;
          const double xi_69 = u_1 * 0.25;
          const double xi_70 = xi_261 * xi_69;
          const double xi_74 = u_1 * xi_71;
          const double xi_75 = xi_261 * xi_74;
          const double xi_76 = xi_68 * -1.0 + xi_70 * -1.0 + xi_73 + xi_75;
          const double xi_78 = xi_68 + xi_70 + xi_73 * -1.0 + xi_75 * -1.0;
          const double xi_86 = xi_259 * xi_69;
          const double xi_88 = xi_259 * xi_74;
          const double xi_93 = xi_45 * -1.0;
          const double xi_112 = rho * u_1 * u_1;
          const double xi_121 = xi_112 + xi_120 + xi_9;
          const double xi_197 = rho * u_1;
          const double xi_199 = xi_198 * (u_0 * xi_197 + xi_120 + xi_247 + xi_252);
          const double xi_200 = xi_196 * -1.0 + xi_199 * -1.0;
          const double xi_201 = xi_196 + xi_199;
          const double u_2 = xi_0 * (vel2Term + xi_21 + xi_24 + xi_260) + xi_259 * xi_7;
          const double xi_27 = u_2 * xi_259;
          const double xi_33 = xi_27 * 0.16666666666666666;
          const double xi_34 = xi_27 * 0.083333333333333329;
          const double xi_35 = omega_shear * xi_34;
          const double xi_36 = xi_33 * -1.0 + xi_35;
          const double xi_41 = omega_shear * xi_32 * -1.0 + xi_26 * 0.33333333333333331 + xi_36 + xi_40;
          const double xi_48 = omega_shear * xi_37 * -1.0 + xi_25 * 0.33333333333333331 + xi_36 + xi_47;
          const double xi_52 = omega_shear * xi_33 * -1.0 + xi_27 * 0.33333333333333331 + xi_40 + xi_47;
          const double xi_58 = xi_34 * -1.0;
          const double xi_63 = xi_27 * xi_60 * -1.0;
          const double xi_64 = xi_26 * xi_55 * -1.0 + xi_32 + xi_61 + xi_62 + xi_63;
          const double xi_65 = xi_30 + xi_59 * -1.0 + xi_64;
          const double xi_66 = xi_35 + xi_58 + xi_65;
          const double xi_79 = xi_30 * -1.0 + xi_59 + xi_64;
          const double xi_80 = xi_35 + xi_58 + xi_79;
          const double xi_82 = xi_27 * xi_55 * -1.0 + xi_33;
          const double xi_83 = xi_50 + xi_81 * -1.0 + xi_82;
          const double xi_85 = xi_39 + xi_65 + xi_84;
          const double xi_87 = u_2 * xi_67;
          const double xi_89 = u_2 * xi_72;
          const double xi_90 = xi_86 + xi_87 + xi_88 * -1.0 + xi_89 * -1.0;
          const double xi_91 = xi_39 + xi_79 + xi_84;
          const double xi_92 = xi_86 * -1.0 + xi_87 * -1.0 + xi_88 + xi_89;
          const double xi_94 = xi_46 + xi_61 + xi_62 + xi_63 + xi_83 + xi_93;
          const double xi_97 = u_2 * xi_261;
          const double xi_98 = xi_97 * 0.25;
          const double xi_100 = xi_71 * xi_97;
          const double xi_101 = xi_100 + xi_96 * -1.0 + xi_98 * -1.0 + xi_99;
          const double xi_102 = xi_100 * -1.0 + xi_96 + xi_98 + xi_99 * -1.0;
          const double xi_103 = xi_50 * -1.0 + xi_81 + xi_82;
          const double xi_104 = xi_103 + xi_46 + xi_61 + xi_62 + xi_63 + xi_93;
          const double xi_115 = rho * u_2 * u_2;
          const double xi_116 = xi_114 + xi_115 * 0.66666666666666663 + xi_245 * 3.0 + xi_265 * 3.0;
          const double xi_117 = omega_even * (xi_112 * 0.66666666666666663 + xi_113 * 1.6666666666666667 + xi_116 + xi_246 * -3.0 + xi_248 * -3.0 + xi_251 * 3.0 + xi_255 * -3.0 + xi_263 * -3.0 + xi_264 * 3.0);
          const double xi_124 = omega_bulk * (xi_113 + xi_115 + xi_119 + xi_121 + xi_123 + xi_17 + xi_22 + xi_253);
          const double xi_127 = omega_even * (xi_112 * 2.3333333333333335 + xi_116 + xi_126 + xi_249 * -5.0 + xi_251 * -2.0 + xi_256 * -5.0 + xi_260 * -5.0 + xi_262 * -5.0 + xi_264 * -2.0);
          const double xi_131 = omega_even * (xi_114 + xi_115 * 3.0 + xi_126 + xi_128 + xi_129 + xi_130 + xi_245 * -4.0 + xi_247 * -7.0 + xi_250 * -7.0 + xi_251 * 5.0 + xi_252 * -7.0 + xi_257 * -7.0 + xi_264 * 5.0 + xi_265 * -4.0);
          const double xi_156 = xi_115 * -1.0 + xi_265;
          const double xi_157 = omega_shear * (xi_121 + xi_156 + xi_16 + xi_2 + xi_251 * -1.0 + xi_256 + xi_6);
          const double xi_158 = xi_157 * 0.125;
          const double xi_159 = xi_107 * -0.11904761904761904 + xi_131 * -0.01984126984126984;
          const double xi_160 = omega_shear * (xi_112 * -1.0 + xi_113 * 2.0 + xi_120 + xi_123 + xi_125 + xi_156 + xi_244 * -2.0 + xi_245 + xi_251 + xi_258 * -2.0 + xi_264 + xi_9);
          const double xi_162 = xi_160 * -0.041666666666666664 + xi_161 * -0.16666666666666666;
          const double xi_163 = xi_108 * -0.10000000000000001 + xi_117 * -0.050000000000000003 + xi_162;
          const double xi_164 = xi_111 * 0.028571428571428571 + xi_127 * 0.014285714285714285 + xi_155 + xi_158 + xi_159 + xi_163;
          const double xi_176 = xi_111 * -0.071428571428571425 + xi_127 * -0.035714285714285712 + xi_159 + xi_160 * 0.083333333333333329 + xi_161 * 0.33333333333333331;
          const double xi_187 = xi_107 * 0.095238095238095233 + xi_111 * -0.042857142857142858 + xi_127 * -0.021428571428571429 + xi_131 * 0.015873015873015872 + xi_155 * -1.0 + xi_158 * -1.0 + xi_163;
          const double xi_190 = xi_157 * 0.0625;
          const double xi_191 = xi_131 * 0.013888888888888888;
          const double xi_193 = xi_110 * 0.083333333333333329 + xi_124 * 0.041666666666666664;
          const double xi_194 = xi_160 * 0.020833333333333332 + xi_161 * 0.083333333333333329 + xi_193;
          const double xi_195 = xi_165 + xi_189 + xi_190 + xi_191 + xi_192 + xi_194;
          const double xi_202 = xi_151 + xi_189 + xi_190 + xi_191 + xi_192 + xi_194;
          const double xi_204 = xi_127 * -0.0071428571428571426;
          const double xi_205 = xi_198 * (u_2 * xi_197 + xi_132 + xi_17 + xi_246);
          const double xi_206 = xi_117 * 0.025000000000000001;
          const double xi_209 = xi_107 * -0.023809523809523808 + xi_131 * -0.003968253968253968;
          const double xi_210 = xi_162 + xi_193 + xi_203 + xi_204 + xi_205 + xi_206 + xi_207 + xi_208 + xi_209;
          const double xi_226 = xi_162 + xi_193 + xi_203 * -1.0 + xi_204 + xi_205 * -1.0 + xi_206 + xi_207 + xi_208 + xi_209;
          const double xi_228 = xi_190 * -1.0;
          const double xi_229 = xi_127 * 0.017857142857142856;
          const double xi_231 = xi_188 + xi_194 + xi_209 + xi_227 + xi_228 + xi_229 + xi_230;
          const double xi_233 = xi_198 * (rho * u_0 * u_2 + xi_10 + xi_166 + xi_249);
          const double xi_234 = xi_232 * -1.0 + xi_233 * -1.0;
          const double xi_240 = xi_232 + xi_233;
          const double xi_243 = xi_186 + xi_194 + xi_209 + xi_227 + xi_228 + xi_229 + xi_230;
          const double forceTerm_0 = xi_25 * xi_28 + xi_25 * -1.0 + xi_26 * xi_28 + xi_26 * -1.0 + xi_27 * xi_28 + xi_27 * -1.0;
          const double forceTerm_1 = xi_29 + xi_31 * -1.0 + xi_41;
          const double forceTerm_2 = xi_29 * -1.0 + xi_31 + xi_41;
          const double forceTerm_3 = xi_42 * -1.0 + xi_44 + xi_48;
          const double forceTerm_4 = xi_42 + xi_44 * -1.0 + xi_48;
          const double forceTerm_5 = xi_49 + xi_51 * -1.0 + xi_52;
          const double forceTerm_6 = xi_49 * -1.0 + xi_51 + xi_52;
          const double forceTerm_7 = xi_57 + xi_66 + xi_76;
          const double forceTerm_8 = xi_66 + xi_77 + xi_78;
          const double forceTerm_9 = xi_57 + xi_78 + xi_80;
          const double forceTerm_10 = xi_76 + xi_77 + xi_80;
          const double forceTerm_11 = xi_83 + xi_85 + xi_90;
          const double forceTerm_12 = xi_83 + xi_91 + xi_92;
          const double forceTerm_13 = xi_101 + xi_57 + xi_94;
          const double forceTerm_14 = xi_102 + xi_77 + xi_94;
          const double forceTerm_15 = xi_103 + xi_85 + xi_92;
          const double forceTerm_16 = xi_103 + xi_90 + xi_91;
          const double forceTerm_17 = xi_102 + xi_104 + xi_57;
          const double forceTerm_18 = xi_101 + xi_104 + xi_77;
          _data_pdfs_20_30_10[ctr_0] = forceTerm_0 + xi_107 * 0.14285714285714285 + xi_108 * 0.20000000000000001 + xi_110 * -1.0 + xi_111 * 0.085714285714285715 + xi_117 * 0.10000000000000001 + xi_124 * -0.5 + xi_127 * 0.042857142857142858 + xi_131 * 0.023809523809523808 + xi_253;
          _data_pdfs_20_31_10[ctr_0] = forceTerm_1 + xi_135 * -1.0 + xi_140 * -1.0 + xi_151 + xi_164 + xi_251;
          _data_pdfs_20_32_10[ctr_0] = forceTerm_2 + xi_135 + xi_140 + xi_164 + xi_165 + xi_264;
          _data_pdfs_20_33_10[ctr_0] = forceTerm_3 + xi_168 + xi_170 + xi_175 + xi_176 + xi_258;
          _data_pdfs_20_34_10[ctr_0] = forceTerm_4 + xi_168 * -1.0 + xi_170 * -1.0 + xi_176 + xi_177 + xi_244;
          _data_pdfs_20_35_10[ctr_0] = forceTerm_5 + xi_179 * -1.0 + xi_181 * -1.0 + xi_186 + xi_187 + xi_265;
          _data_pdfs_20_36_10[ctr_0] = forceTerm_6 + xi_179 + xi_181 + xi_187 + xi_188 + xi_245;
          _data_pdfs_20_37_10[ctr_0] = forceTerm_7 + xi_177 + xi_195 + xi_200 + xi_252;
          _data_pdfs_20_38_10[ctr_0] = forceTerm_8 + xi_175 + xi_195 + xi_201 + xi_257;
          _data_pdfs_20_39_10[ctr_0] = forceTerm_9 + xi_177 + xi_201 + xi_202 + xi_250;
          _data_pdfs_20_310_10[ctr_0] = forceTerm_10 + xi_175 + xi_200 + xi_202 + xi_247;
          _data_pdfs_20_311_10[ctr_0] = forceTerm_11 + xi_210 + xi_219 + xi_224 + xi_263;
          _data_pdfs_20_312_10[ctr_0] = forceTerm_12 + xi_219 + xi_225 + xi_226 + xi_248;
          _data_pdfs_20_313_10[ctr_0] = forceTerm_13 + xi_231 + xi_234 + xi_239 + xi_256;
          _data_pdfs_20_314_10[ctr_0] = forceTerm_14 + xi_231 + xi_240 + xi_241 + xi_260;
          _data_pdfs_20_315_10[ctr_0] = forceTerm_15 + xi_224 + xi_226 + xi_242 + xi_246;
          _data_pdfs_20_316_10[ctr_0] = forceTerm_16 + xi_210 + xi_225 + xi_242 + xi_255;
          _data_pdfs_20_317_10[ctr_0] = forceTerm_17 + xi_239 + xi_240 + xi_243 + xi_262;
          _data_pdfs_20_318_10[ctr_0] = forceTerm_18 + xi_234 + xi_241 + xi_243 + xi_249;
        }
      }
    }
  }
}
} // namespace internal_25bc51f30ec2c20f3ee9796f7dcb65c6

void CollideSweepDoublePrecisionThermalizedAVX::run(IBlock *block) {
  auto force = block->getData<field::GhostLayerField<double, 3>>(forceID);
  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);

  auto &omega_bulk = this->omega_bulk_;
  auto block_offset_1 = this->block_offset_1_;
  auto &seed = this->seed_;
  auto &omega_even = this->omega_even_;
  auto &kT = this->kT_;
  auto &omega_odd = this->omega_odd_;
  auto block_offset_2 = this->block_offset_2_;
  auto &time_step = this->time_step_;
  auto block_offset_0 = this->block_offset_0_;
  auto &omega_shear = this->omega_shear_;
  block_offset_generator(block, block_offset_0, block_offset_1, block_offset_2);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force->nrOfGhostLayers()));
  double *RESTRICT const _data_force = force->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()));
  double *RESTRICT _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
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
  internal_25bc51f30ec2c20f3ee9796f7dcb65c6::collidesweepdoubleprecisionthermalizedavx_collidesweepdoubleprecisionthermalizedavx(_data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, block_offset_0, block_offset_1, block_offset_2, kT, omega_bulk, omega_even, omega_odd, omega_shear, seed, time_step);
}

void CollideSweepDoublePrecisionThermalizedAVX::runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks, const CellInterval &globalCellInterval, cell_idx_t ghostLayers, IBlock *block) {
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
  auto block_offset_1 = this->block_offset_1_;
  auto &seed = this->seed_;
  auto &omega_even = this->omega_even_;
  auto &kT = this->kT_;
  auto &omega_odd = this->omega_odd_;
  auto block_offset_2 = this->block_offset_2_;
  auto &time_step = this->time_step_;
  auto block_offset_0 = this->block_offset_0_;
  auto &omega_shear = this->omega_shear_;
  block_offset_generator(block, block_offset_0, block_offset_1, block_offset_2);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(force->nrOfGhostLayers()));
  double *RESTRICT const _data_force = force->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
  double *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
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
  internal_25bc51f30ec2c20f3ee9796f7dcb65c6::collidesweepdoubleprecisionthermalizedavx_collidesweepdoubleprecisionthermalizedavx(_data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, block_offset_0, block_offset_1, block_offset_2, kT, omega_bulk, omega_even, omega_odd, omega_shear, seed, time_step);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif