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
//! \\file StreamCollideSweepThermalizedSinglePrecisionAVX.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 007e77e077ad9d22b5eed6f3d3118240993e553c

#include <cmath>

#include "StreamCollideSweepThermalizedSinglePrecisionAVX.h"
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

namespace internal_482af0085ade5d86c474548ac4b32084 {
static FUNC_PREFIX void streamcollidesweepthermalizedsingleprecisionavx_streamcollidesweepthermalizedsingleprecisionavx(float *RESTRICT const _data_force, float *RESTRICT const _data_pdfs, float *RESTRICT _data_pdfs_tmp, int64_t const _size_force_0, int64_t const _size_force_1, int64_t const _size_force_2, int64_t const _stride_force_1, int64_t const _stride_force_2, int64_t const _stride_force_3, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, int64_t const _stride_pdfs_tmp_1, int64_t const _stride_pdfs_tmp_2, int64_t const _stride_pdfs_tmp_3, uint32_t block_offset_0, uint32_t block_offset_1, uint32_t block_offset_2, float kT, float omega_bulk, float omega_even, float omega_odd, float omega_shear, uint32_t seed, uint32_t time_step) {
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    const float xi_20 = omega_bulk * 0.5f;
    const float xi_47 = omega_shear * 0.041666666666666664f;
    const float xi_51 = omega_bulk * 0.041666666666666664f;
    const float xi_62 = omega_shear * 0.125f;
    const float xi_97 = 3.7416573867739413f;
    const float xi_100 = 5.4772255750516612f;
    const float xi_104 = 2.4494897427831779f;
    const float xi_107 = 8.3666002653407556f;
    const float xi_150 = omega_odd * 0.25f;
    const float xi_160 = omega_odd * 0.083333333333333329f;
    const float xi_173 = 1.7320508075688772f;
    const float xi_217 = omega_shear * 0.25f;
    const float xi_223 = omega_odd * 0.041666666666666664f;
    const float xi_226 = omega_odd * 0.125f;
    const float rr_0 = 0.0f;
    const float xi_45 = rr_0 * 0.041666666666666664f;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int64_t ctr_2 = 1; ctr_2 < _size_force_2 - 1; ctr_2 += 1) {
      for (int64_t ctr_1 = 1; ctr_1 < _size_force_1 - 1; ctr_1 += 1) {
        {
          for (int64_t ctr_0 = 1; ctr_0 < (int64_t)((_size_force_0 - 2) / (8)) * (8) + 1; ctr_0 += 8) {

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
            const __m256 xi_2 = _mm256_add_ps(_mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3 + ctr_0 - 1]), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3 + ctr_0 - 1]));
            const __m256 xi_3 = _mm256_add_ps(_mm256_add_ps(xi_2, _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3 + ctr_0 - 1])), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1]));
            const __m256 xi_4 = _mm256_add_ps(_mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3 + ctr_0]), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3 + ctr_0 + 1]));
            const __m256 xi_5 = _mm256_add_ps(_mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0]), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0]));
            const __m256 xi_6 = _mm256_add_ps(_mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3 + ctr_0]), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0]));
            const __m256 xi_7 = _mm256_add_ps(_mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3 + ctr_0 + 1]), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1]));
            const __m256 xi_8 = _mm256_add_ps(xi_7, _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1]));
            const __m256 xi_9 = _mm256_add_ps(_mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0]), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0]));
            const __m256 xi_11 = _mm256_add_ps(_mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3 + ctr_0 + 1]), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1]));
            const __m256 xi_12 = _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3 + ctr_0 - 1])), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3 + ctr_0 - 1]));
            const __m256 xi_13 = _mm256_add_ps(_mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1]), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0]));
            const __m256 xi_14 = _mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0]));
            const __m256 xi_15 = _mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3 + ctr_0 - 1]));
            const __m256 xi_16 = _mm256_add_ps(_mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1]), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0]));
            const __m256 xi_21 = _mm256_mul_ps(_mm256_set_ps(0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f), _mm256_loadu_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0]));
            const __m256 xi_22 = _mm256_mul_ps(_mm256_set_ps(0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f), _mm256_loadu_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0]));
            const __m256 xi_33 = _mm256_mul_ps(_mm256_set_ps(0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f), _mm256_load_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0]));
            const __m256 xi_34 = _mm256_mul_ps(_mm256_set_ps(0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f), _mm256_load_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0]));
            const __m256 xi_39 = _mm256_mul_ps(_mm256_set_ps(0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f), _mm256_loadu_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0]));
            const __m256 xi_40 = _mm256_mul_ps(_mm256_set_ps(0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f), _mm256_loadu_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0]));
            const __m256 xi_58 = _mm256_mul_ps(_mm256_set_ps(0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f), _mm256_loadu_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0]));
            const __m256 xi_63 = _mm256_mul_ps(_mm256_set_ps(xi_62, xi_62, xi_62, xi_62, xi_62, xi_62, xi_62, xi_62), _mm256_loadu_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0]));
            const __m256 xi_98 = _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_3_0);
            const __m256 xi_101 = _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_3_2);
            const __m256 xi_103 = _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_1_1);
            const __m256 xi_108 = _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_3_1);
            const __m256 xi_112 = _mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_load_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + ctr_0]));
            const __m256 xi_118 = _mm256_add_ps(xi_11, xi_3);
            const __m256 xi_121 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0])), _mm256_mul_ps(_mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_ps(_mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_ps(_mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_ps(_mm256_set_ps(5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3 + ctr_0 + 1]))), _mm256_mul_ps(_mm256_set_ps(5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3 + ctr_0 - 1])));
            const __m256 xi_124 = _mm256_mul_ps(_mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1]));
            const __m256 xi_125 = _mm256_mul_ps(_mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3 + ctr_0 - 1]));
            const __m256 xi_126 = _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1])), _mm256_mul_ps(_mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1])));
            const __m256 xi_129 = _mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3 + ctr_0]));
            const __m256 xi_133 = _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_0_1);
            const __m256 xi_137 = _mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1]));
            const __m256 xi_138 = _mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1]));
            const __m256 xi_139 = _mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3 + ctr_0]));
            const __m256 xi_140 = _mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0]));
            const __m256 xi_145 = _mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0]));
            const __m256 xi_146 = _mm256_add_ps(xi_14, _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0]));
            const __m256 xi_147 = _mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0]));
            const __m256 xi_148 = _mm256_add_ps(xi_147, _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3 + ctr_0]));
            const __m256 xi_149 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_145, xi_146), xi_148), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0]));
            const __m256 xi_151 = _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_2_3);
            const __m256 xi_155 = _mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0]));
            const __m256 xi_156 = _mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0]));
            const __m256 xi_157 = _mm256_mul_ps(_mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3 + ctr_0 - 1]));
            const __m256 xi_158 = _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3 + ctr_0 - 1])), _mm256_mul_ps(_mm256_set_ps(-2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3 + ctr_0 + 1])));
            const __m256 xi_159 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1])), _mm256_mul_ps(xi_157, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_148), xi_155), xi_156), xi_158), xi_5);
            const __m256 xi_161 = _mm256_mul_ps(xi_159, _mm256_set_ps(xi_160, xi_160, xi_160, xi_160, xi_160, xi_160, xi_160, xi_160));
            const __m256 xi_162 = _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_1_2);
            const __m256 xi_174 = _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_0_0);
            const __m256 xi_185 = _mm256_add_ps(xi_15, _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1]));
            const __m256 xi_186 = _mm256_add_ps(xi_138, xi_185);
            const __m256 xi_187 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3 + ctr_0 + 1])), xi_186), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3 + ctr_0 - 1])), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1]));
            const __m256 xi_188 = _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_2_1);
            const __m256 xi_189 = _mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3 + ctr_0 - 1]));
            const __m256 xi_190 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1])), _mm256_mul_ps(xi_157, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_158, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_186, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_189, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_7, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)));
            const __m256 xi_191 = _mm256_mul_ps(xi_190, _mm256_set_ps(xi_160, xi_160, xi_160, xi_160, xi_160, xi_160, xi_160, xi_160));
            const __m256 xi_192 = _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_2_0);
            const __m256 xi_198 = _mm256_add_ps(_mm256_add_ps(xi_145, xi_155), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0]));
            const __m256 xi_199 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_139, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_198, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0])));
            const __m256 xi_200 = _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_2_2);
            const __m256 xi_201 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_124, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_125, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_126), xi_140), xi_198), xi_6);
            const __m256 xi_202 = _mm256_mul_ps(xi_201, _mm256_set_ps(xi_160, xi_160, xi_160, xi_160, xi_160, xi_160, xi_160, xi_160));
            const __m256 xi_203 = _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_1_3);
            const __m256 xi_224 = _mm256_mul_ps(xi_159, _mm256_set_ps(xi_223, xi_223, xi_223, xi_223, xi_223, xi_223, xi_223, xi_223));
            const __m256 xi_227 = _mm256_mul_ps(xi_149, _mm256_set_ps(xi_226, xi_226, xi_226, xi_226, xi_226, xi_226, xi_226, xi_226));
            const __m256 xi_233 = _mm256_mul_ps(xi_201, _mm256_set_ps(xi_223, xi_223, xi_223, xi_223, xi_223, xi_223, xi_223, xi_223));
            const __m256 xi_234 = _mm256_mul_ps(xi_199, _mm256_set_ps(xi_226, xi_226, xi_226, xi_226, xi_226, xi_226, xi_226, xi_226));
            const __m256 xi_254 = _mm256_mul_ps(xi_187, _mm256_set_ps(xi_226, xi_226, xi_226, xi_226, xi_226, xi_226, xi_226, xi_226));
            const __m256 xi_255 = _mm256_mul_ps(xi_190, _mm256_set_ps(xi_223, xi_223, xi_223, xi_223, xi_223, xi_223, xi_223, xi_223));
            const __m256 xi_23 = _mm256_mul_ps(xi_22, _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0, rr_0, rr_0, rr_0));
            const __m256 xi_35 = _mm256_mul_ps(xi_34, _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0, rr_0, rr_0, rr_0));
            const __m256 xi_41 = _mm256_mul_ps(xi_40, _mm256_set_ps(rr_0, rr_0, rr_0, rr_0, rr_0, rr_0, rr_0, rr_0));
            const __m256 xi_46 = _mm256_mul_ps(_mm256_set_ps(xi_45, xi_45, xi_45, xi_45, xi_45, xi_45, xi_45, xi_45), _mm256_load_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0]));
            const __m256 xi_50 = _mm256_mul_ps(_mm256_set_ps(xi_45, xi_45, xi_45, xi_45, xi_45, xi_45, xi_45, xi_45), _mm256_loadu_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0]));
            const __m256 xi_72 = _mm256_mul_ps(_mm256_set_ps(xi_45, xi_45, xi_45, xi_45, xi_45, xi_45, xi_45, xi_45), _mm256_loadu_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0]));
            const __m256 vel0Term = _mm256_add_ps(xi_3, _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3 + ctr_0 - 1]));
            const __m256 vel1Term = _mm256_add_ps(xi_4, xi_5);
            const __m256 vel2Term = _mm256_add_ps(xi_6, _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1]));
            const __m256 delta_rho = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(vel0Term, vel1Term), vel2Term), xi_8), xi_9), _mm256_load_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + ctr_0])), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0]));
            const __m256 rho = _mm256_add_ps(_mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f), delta_rho);
            const __m256 xi_95 = _mm256_mul_ps(rho, _mm256_set_ps(kT, kT, kT, kT, kT, kT, kT, kT));
            const __m256 xi_96 = _mm256_sqrt_ps(_mm256_mul_ps(xi_95, _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_mul_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_set_ps(omega_even, omega_even, omega_even, omega_even, omega_even, omega_even, omega_even, omega_even)), _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f)), _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_set_ps(omega_even, omega_even, omega_even, omega_even, omega_even, omega_even, omega_even, omega_even)), _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f)))), _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f))));
            const __m256 xi_99 = _mm256_mul_ps(_mm256_mul_ps(xi_96, xi_98), _mm256_set_ps(xi_97, xi_97, xi_97, xi_97, xi_97, xi_97, xi_97, xi_97));
            const __m256 xi_102 = _mm256_mul_ps(_mm256_mul_ps(xi_101, xi_96), _mm256_set_ps(xi_100, xi_100, xi_100, xi_100, xi_100, xi_100, xi_100, xi_100));
            const __m256 xi_105 = _mm256_sqrt_ps(_mm256_mul_ps(xi_95, _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_mul_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk)), _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f)), _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk)), _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f)))), _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f))));
            const __m256 xi_106 = _mm256_mul_ps(_mm256_mul_ps(xi_103, xi_105), _mm256_set_ps(xi_104, xi_104, xi_104, xi_104, xi_104, xi_104, xi_104, xi_104));
            const __m256 xi_109 = _mm256_mul_ps(_mm256_mul_ps(xi_108, xi_96), _mm256_set_ps(xi_107, xi_107, xi_107, xi_107, xi_107, xi_107, xi_107, xi_107));
            const __m256 xi_131 = _mm256_mul_ps(xi_99, _mm256_set_ps(0.11904761904761904f, 0.11904761904761904f, 0.11904761904761904f, 0.11904761904761904f, 0.11904761904761904f, 0.11904761904761904f, 0.11904761904761904f, 0.11904761904761904f));
            const __m256 xi_134 = _mm256_sqrt_ps(_mm256_mul_ps(xi_95, _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_mul_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear)), _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f)), _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear)), _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f)))), _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f))));
            const __m256 xi_135 = _mm256_mul_ps(xi_134, _mm256_set_ps(0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f));
            const __m256 xi_136 = _mm256_mul_ps(xi_133, xi_135);
            const __m256 xi_152 = _mm256_sqrt_ps(_mm256_mul_ps(xi_95, _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_mul_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_set_ps(omega_odd, omega_odd, omega_odd, omega_odd, omega_odd, omega_odd, omega_odd, omega_odd)), _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f)), _mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_set_ps(omega_odd, omega_odd, omega_odd, omega_odd, omega_odd, omega_odd, omega_odd, omega_odd)), _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f)))), _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f))));
            const __m256 xi_153 = _mm256_mul_ps(xi_152, _mm256_set_ps(1.4142135623730951f, 1.4142135623730951f, 1.4142135623730951f, 1.4142135623730951f, 1.4142135623730951f, 1.4142135623730951f, 1.4142135623730951f, 1.4142135623730951f));
            const __m256 xi_154 = _mm256_mul_ps(xi_153, _mm256_set_ps(0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f));
            const __m256 xi_163 = _mm256_mul_ps(xi_152, _mm256_set_ps(xi_104, xi_104, xi_104, xi_104, xi_104, xi_104, xi_104, xi_104));
            const __m256 xi_164 = _mm256_mul_ps(xi_163, _mm256_set_ps(0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f));
            const __m256 xi_165 = _mm256_mul_ps(xi_162, xi_164);
            const __m256 xi_166 = _mm256_add_ps(xi_161, xi_165);
            const __m256 xi_167 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_149, _mm256_set_ps(xi_150, xi_150, xi_150, xi_150, xi_150, xi_150, xi_150, xi_150)), _mm256_mul_ps(xi_151, xi_154)), xi_166);
            const __m256 xi_169 = _mm256_mul_ps(xi_102, _mm256_set_ps(0.10000000000000001f, 0.10000000000000001f, 0.10000000000000001f, 0.10000000000000001f, 0.10000000000000001f, 0.10000000000000001f, 0.10000000000000001f, 0.10000000000000001f));
            const __m256 xi_175 = _mm256_mul_ps(_mm256_mul_ps(xi_134, xi_174), _mm256_set_ps(xi_173, xi_173, xi_173, xi_173, xi_173, xi_173, xi_173, xi_173));
            const __m256 xi_176 = _mm256_mul_ps(xi_175, _mm256_set_ps(0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f));
            const __m256 xi_184 = _mm256_mul_ps(xi_109, _mm256_set_ps(0.071428571428571425f, 0.071428571428571425f, 0.071428571428571425f, 0.071428571428571425f, 0.071428571428571425f, 0.071428571428571425f, 0.071428571428571425f, 0.071428571428571425f));
            const __m256 xi_193 = _mm256_mul_ps(xi_164, xi_192);
            const __m256 xi_194 = _mm256_add_ps(xi_191, xi_193);
            const __m256 xi_195 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_187, _mm256_set_ps(xi_150, xi_150, xi_150, xi_150, xi_150, xi_150, xi_150, xi_150)), _mm256_mul_ps(xi_154, xi_188)), xi_194);
            const __m256 xi_197 = _mm256_mul_ps(xi_109, _mm256_set_ps(0.042857142857142858f, 0.042857142857142858f, 0.042857142857142858f, 0.042857142857142858f, 0.042857142857142858f, 0.042857142857142858f, 0.042857142857142858f, 0.042857142857142858f));
            const __m256 xi_204 = _mm256_mul_ps(xi_164, xi_203);
            const __m256 xi_205 = _mm256_add_ps(xi_202, xi_204);
            const __m256 xi_206 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_199, _mm256_set_ps(xi_150, xi_150, xi_150, xi_150, xi_150, xi_150, xi_150, xi_150)), _mm256_mul_ps(xi_154, xi_200)), xi_205);
            const __m256 xi_207 = _mm256_mul_ps(_mm256_mul_ps(xi_133, xi_134), _mm256_set_ps(0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f));
            const __m256 xi_210 = _mm256_mul_ps(xi_99, _mm256_set_ps(0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f));
            const __m256 xi_214 = _mm256_add_ps(_mm256_mul_ps(xi_191, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_193, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)));
            const __m256 xi_215 = _mm256_mul_ps(xi_135, _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_0_2));
            const __m256 xi_222 = _mm256_mul_ps(xi_135, _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_1_0));
            const __m256 xi_228 = _mm256_mul_ps(xi_163, _mm256_set_ps(0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f));
            const __m256 xi_229 = _mm256_mul_ps(xi_162, xi_228);
            const __m256 xi_230 = _mm256_mul_ps(xi_153, _mm256_set_ps(0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f));
            const __m256 xi_231 = _mm256_mul_ps(xi_151, xi_230);
            const __m256 xi_235 = _mm256_mul_ps(xi_203, xi_228);
            const __m256 xi_236 = _mm256_mul_ps(xi_200, xi_230);
            const __m256 xi_237 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_233, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_235, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_234), xi_236);
            const __m256 xi_239 = _mm256_mul_ps(xi_109, _mm256_set_ps(0.014285714285714285f, 0.014285714285714285f, 0.014285714285714285f, 0.014285714285714285f, 0.014285714285714285f, 0.014285714285714285f, 0.014285714285714285f, 0.014285714285714285f));
            const __m256 xi_241 = _mm256_mul_ps(xi_99, _mm256_set_ps(0.023809523809523808f, 0.023809523809523808f, 0.023809523809523808f, 0.023809523809523808f, 0.023809523809523808f, 0.023809523809523808f, 0.023809523809523808f, 0.023809523809523808f));
            const __m256 xi_244 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_234, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_236, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_233), xi_235);
            const __m256 xi_246 = _mm256_mul_ps(xi_207, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f));
            const __m256 xi_249 = _mm256_mul_ps(xi_109, _mm256_set_ps(0.035714285714285712f, 0.035714285714285712f, 0.035714285714285712f, 0.035714285714285712f, 0.035714285714285712f, 0.035714285714285712f, 0.035714285714285712f, 0.035714285714285712f));
            const __m256 xi_251 = _mm256_mul_ps(xi_135, _mm256_add_ps(_mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f), random_0_3));
            const __m256 xi_256 = _mm256_mul_ps(xi_188, xi_230);
            const __m256 xi_257 = _mm256_mul_ps(xi_192, xi_228);
            const __m256 xi_258 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_254, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_256, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_255), xi_257);
            const __m256 xi_260 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_255, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_257, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_254), xi_256);
            const __m256 xi_0 = _mm256_div_ps(_mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f), rho);
            const __m256 xi_10 = _mm256_mul_ps(xi_0, _mm256_set_ps(0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f));
            const __m256 u_0 = _mm256_add_ps(_mm256_mul_ps(xi_0, _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_11, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_8, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), vel0Term)), _mm256_mul_ps(xi_10, _mm256_load_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0])));
            const __m256 xi_17 = _mm256_mul_ps(u_0, _mm256_load_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0]));
            const __m256 xi_28 = _mm256_mul_ps(xi_17, _mm256_set_ps(0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f));
            const __m256 xi_29 = _mm256_mul_ps(xi_28, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f));
            const __m256 xi_30 = _mm256_mul_ps(xi_17, _mm256_set_ps(0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f));
            const __m256 xi_31 = _mm256_add_ps(_mm256_mul_ps(xi_30, _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear)), xi_29);
            const __m256 xi_48 = _mm256_add_ps(_mm256_mul_ps(xi_17, _mm256_set_ps(xi_47, xi_47, xi_47, xi_47, xi_47, xi_47, xi_47, xi_47)), xi_29);
            const __m256 xi_49 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_46, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_34), xi_48);
            const __m256 xi_52 = _mm256_mul_ps(xi_17, _mm256_set_ps(xi_51, xi_51, xi_51, xi_51, xi_51, xi_51, xi_51, xi_51));
            const __m256 xi_59 = _mm256_mul_ps(u_0, xi_58);
            const __m256 xi_64 = _mm256_mul_ps(u_0, xi_63);
            const __m256 xi_68 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_34, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_46), xi_48);
            const __m256 xi_75 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(u_0, _mm256_set_ps(-0.083333333333333329f, -0.083333333333333329f, -0.083333333333333329f, -0.083333333333333329f, -0.083333333333333329f, -0.083333333333333329f, -0.083333333333333329f, -0.083333333333333329f)), _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear)), _mm256_load_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0]));
            const __m256 xi_85 = _mm256_mul_ps(u_0, _mm256_loadu_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0]));
            const __m256 xi_86 = _mm256_mul_ps(xi_85, _mm256_set_ps(0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f));
            const __m256 xi_89 = _mm256_mul_ps(xi_85, _mm256_set_ps(xi_62, xi_62, xi_62, xi_62, xi_62, xi_62, xi_62, xi_62));
            const __m256 xi_111 = _mm256_mul_ps(u_0, u_0);
            const __m256 u_1 = _mm256_add_ps(_mm256_mul_ps(xi_0, _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_12, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_13, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_9, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), vel1Term)), _mm256_mul_ps(xi_10, _mm256_loadu_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0])));
            const __m256 xi_18 = _mm256_mul_ps(u_1, _mm256_loadu_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0]));
            const __m256 xi_26 = _mm256_mul_ps(xi_18, _mm256_set_ps(0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f));
            const __m256 xi_36 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(u_1, _mm256_set_ps(-0.083333333333333329f, -0.083333333333333329f, -0.083333333333333329f, -0.083333333333333329f, -0.083333333333333329f, -0.083333333333333329f, -0.083333333333333329f, -0.083333333333333329f)), _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear)), _mm256_loadu_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0]));
            const __m256 xi_42 = _mm256_mul_ps(xi_26, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f));
            const __m256 xi_43 = _mm256_mul_ps(xi_18, _mm256_set_ps(0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f));
            const __m256 xi_53 = _mm256_mul_ps(xi_18, _mm256_set_ps(xi_51, xi_51, xi_51, xi_51, xi_51, xi_51, xi_51, xi_51));
            const __m256 xi_60 = _mm256_mul_ps(u_1, _mm256_set_ps(0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f));
            const __m256 xi_61 = _mm256_mul_ps(xi_60, _mm256_load_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0]));
            const __m256 xi_65 = _mm256_mul_ps(u_1, _mm256_set_ps(xi_62, xi_62, xi_62, xi_62, xi_62, xi_62, xi_62, xi_62));
            const __m256 xi_66 = _mm256_mul_ps(xi_65, _mm256_load_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0]));
            const __m256 xi_67 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_64, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_66, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_59), xi_61);
            const __m256 xi_69 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_59, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_61, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_64), xi_66);
            const __m256 xi_77 = _mm256_mul_ps(xi_60, _mm256_loadu_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0]));
            const __m256 xi_79 = _mm256_mul_ps(xi_65, _mm256_loadu_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0]));
            const __m256 xi_110 = _mm256_mul_ps(rho, _mm256_mul_ps(u_1, u_1));
            const __m256 xi_117 = _mm256_mul_ps(xi_110, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f));
            const __m256 xi_216 = _mm256_mul_ps(rho, u_1);
            const __m256 xi_218 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1])), _mm256_mul_ps(u_0, xi_216)), xi_12), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3 + ctr_0 + 1])), _mm256_set_ps(xi_217, xi_217, xi_217, xi_217, xi_217, xi_217, xi_217, xi_217));
            const __m256 xi_219 = _mm256_add_ps(_mm256_mul_ps(xi_215, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_218, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)));
            const __m256 xi_220 = _mm256_add_ps(xi_215, xi_218);
            const __m256 u_2 = _mm256_add_ps(_mm256_mul_ps(xi_0, _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_14, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_15, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_16, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1]))), _mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0]))), vel2Term)), _mm256_mul_ps(xi_10, _mm256_loadu_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0])));
            const __m256 xi_19 = _mm256_mul_ps(u_2, _mm256_loadu_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0]));
            const __m256 xi_24 = _mm256_mul_ps(xi_19, _mm256_set_ps(0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f, 0.16666666666666666f));
            const __m256 xi_25 = _mm256_mul_ps(xi_24, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f));
            const __m256 xi_27 = _mm256_mul_ps(xi_19, _mm256_set_ps(0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f));
            const __m256 xi_32 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_18, _mm256_set_ps(0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f)), _mm256_mul_ps(xi_27, _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear))), _mm256_mul_ps(_mm256_mul_ps(xi_26, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear))), xi_25), xi_31);
            const __m256 xi_37 = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(u_2, _mm256_set_ps(-0.083333333333333329f, -0.083333333333333329f, -0.083333333333333329f, -0.083333333333333329f, -0.083333333333333329f, -0.083333333333333329f, -0.083333333333333329f, -0.083333333333333329f)), _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear)), _mm256_loadu_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0]));
            const __m256 xi_38 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_28, _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear)), _mm256_mul_ps(_mm256_mul_ps(u_0, _mm256_set_ps(-0.33333333333333331f, -0.33333333333333331f, -0.33333333333333331f, -0.33333333333333331f, -0.33333333333333331f, -0.33333333333333331f, -0.33333333333333331f, -0.33333333333333331f)), _mm256_load_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0]))), xi_24), xi_26), xi_36), xi_37);
            const __m256 xi_44 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_19, _mm256_set_ps(0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f)), _mm256_mul_ps(xi_43, _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear))), _mm256_mul_ps(_mm256_mul_ps(xi_24, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear))), xi_31), xi_42);
            const __m256 xi_54 = _mm256_mul_ps(xi_19, _mm256_set_ps(xi_51, xi_51, xi_51, xi_51, xi_51, xi_51, xi_51, xi_51));
            const __m256 xi_55 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_18, _mm256_set_ps(xi_47, xi_47, xi_47, xi_47, xi_47, xi_47, xi_47, xi_47)), xi_42), xi_52), xi_53), xi_54);
            const __m256 xi_56 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_22, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_50), xi_55);
            const __m256 xi_57 = _mm256_add_ps(_mm256_add_ps(xi_27, xi_37), xi_56);
            const __m256 xi_70 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_50, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_22), xi_55);
            const __m256 xi_71 = _mm256_add_ps(_mm256_add_ps(xi_27, xi_37), xi_70);
            const __m256 xi_73 = _mm256_add_ps(_mm256_mul_ps(xi_19, _mm256_set_ps(xi_47, xi_47, xi_47, xi_47, xi_47, xi_47, xi_47, xi_47)), xi_25);
            const __m256 xi_74 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_40, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_72), xi_73);
            const __m256 xi_76 = _mm256_add_ps(_mm256_add_ps(xi_30, xi_56), xi_75);
            const __m256 xi_78 = _mm256_mul_ps(u_2, xi_58);
            const __m256 xi_80 = _mm256_mul_ps(u_2, xi_63);
            const __m256 xi_81 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_77, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_78, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_79), xi_80);
            const __m256 xi_82 = _mm256_add_ps(_mm256_add_ps(xi_30, xi_70), xi_75);
            const __m256 xi_83 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_79, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_80, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_77), xi_78);
            const __m256 xi_84 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_36, xi_43), xi_52), xi_53), xi_54), xi_74);
            const __m256 xi_87 = _mm256_mul_ps(u_2, _mm256_load_ps(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0]));
            const __m256 xi_88 = _mm256_mul_ps(xi_87, _mm256_set_ps(0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f, 0.25f));
            const __m256 xi_90 = _mm256_mul_ps(xi_87, _mm256_set_ps(xi_62, xi_62, xi_62, xi_62, xi_62, xi_62, xi_62, xi_62));
            const __m256 xi_91 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_89, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_90, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_86), xi_88);
            const __m256 xi_92 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_86, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_88, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_89), xi_90);
            const __m256 xi_93 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_72, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_40), xi_73);
            const __m256 xi_94 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_36, xi_43), xi_52), xi_53), xi_54), xi_93);
            const __m256 xi_113 = _mm256_mul_ps(rho, _mm256_mul_ps(u_2, u_2));
            const __m256 xi_114 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3 + ctr_0])), _mm256_mul_ps(_mm256_set_ps(3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_ps(xi_113, _mm256_set_ps(0.66666666666666663f, 0.66666666666666663f, 0.66666666666666663f, 0.66666666666666663f, 0.66666666666666663f, 0.66666666666666663f, 0.66666666666666663f, 0.66666666666666663f))), xi_112);
            const __m256 xi_115 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3 + ctr_0])), _mm256_mul_ps(_mm256_set_ps(3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_ps(xi_110, _mm256_set_ps(0.66666666666666663f, 0.66666666666666663f, 0.66666666666666663f, 0.66666666666666663f, 0.66666666666666663f, 0.66666666666666663f, 0.66666666666666663f, 0.66666666666666663f))), _mm256_mul_ps(_mm256_set_ps(-3.0f, -3.0f, -3.0f, -3.0f, -3.0f, -3.0f, -3.0f, -3.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_ps(_mm256_set_ps(-3.0f, -3.0f, -3.0f, -3.0f, -3.0f, -3.0f, -3.0f, -3.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_ps(_mm256_set_ps(-3.0f, -3.0f, -3.0f, -3.0f, -3.0f, -3.0f, -3.0f, -3.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_ps(_mm256_set_ps(-3.0f, -3.0f, -3.0f, -3.0f, -3.0f, -3.0f, -3.0f, -3.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_ps(_mm256_mul_ps(rho, xi_111), _mm256_set_ps(1.6666666666666667f, 1.6666666666666667f, 1.6666666666666667f, 1.6666666666666667f, 1.6666666666666667f, 1.6666666666666667f, 1.6666666666666667f, 1.6666666666666667f))), xi_114);
            const __m256 xi_116 = _mm256_mul_ps(xi_115, _mm256_set_ps(omega_even, omega_even, omega_even, omega_even, omega_even, omega_even, omega_even, omega_even));
            const __m256 xi_119 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_112, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_117, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_118, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_13, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_16, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_5, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(rho, xi_111)), xi_113);
            const __m256 xi_120 = _mm256_mul_ps(xi_119, _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk));
            const __m256 xi_122 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_110, _mm256_set_ps(2.3333333333333335f, 2.3333333333333335f, 2.3333333333333335f, 2.3333333333333335f, 2.3333333333333335f, 2.3333333333333335f, 2.3333333333333335f, 2.3333333333333335f)), _mm256_mul_ps(_mm256_set_ps(-2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3 + ctr_0]))), _mm256_mul_ps(_mm256_set_ps(-2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_ps(_mm256_set_ps(-5.0f, -5.0f, -5.0f, -5.0f, -5.0f, -5.0f, -5.0f, -5.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1]))), _mm256_mul_ps(_mm256_set_ps(-5.0f, -5.0f, -5.0f, -5.0f, -5.0f, -5.0f, -5.0f, -5.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1]))), _mm256_mul_ps(_mm256_set_ps(-5.0f, -5.0f, -5.0f, -5.0f, -5.0f, -5.0f, -5.0f, -5.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3 + ctr_0 - 1]))), _mm256_mul_ps(_mm256_set_ps(-5.0f, -5.0f, -5.0f, -5.0f, -5.0f, -5.0f, -5.0f, -5.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1]))), xi_114), xi_121);
            const __m256 xi_123 = _mm256_mul_ps(xi_122, _mm256_set_ps(omega_even, omega_even, omega_even, omega_even, omega_even, omega_even, omega_even, omega_even));
            const __m256 xi_127 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_113, _mm256_set_ps(3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f, 3.0f)), _mm256_mul_ps(_mm256_set_ps(5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3 + ctr_0]))), _mm256_mul_ps(_mm256_set_ps(5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f, 5.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_ps(_mm256_set_ps(-4.0f, -4.0f, -4.0f, -4.0f, -4.0f, -4.0f, -4.0f, -4.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_ps(_mm256_set_ps(-4.0f, -4.0f, -4.0f, -4.0f, -4.0f, -4.0f, -4.0f, -4.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_ps(_mm256_set_ps(-7.0f, -7.0f, -7.0f, -7.0f, -7.0f, -7.0f, -7.0f, -7.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3 + ctr_0 + 1]))), _mm256_mul_ps(_mm256_set_ps(-7.0f, -7.0f, -7.0f, -7.0f, -7.0f, -7.0f, -7.0f, -7.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1]))), _mm256_mul_ps(_mm256_set_ps(-7.0f, -7.0f, -7.0f, -7.0f, -7.0f, -7.0f, -7.0f, -7.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3 + ctr_0 - 1]))), _mm256_mul_ps(_mm256_set_ps(-7.0f, -7.0f, -7.0f, -7.0f, -7.0f, -7.0f, -7.0f, -7.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3 + ctr_0 - 1]))), xi_112), xi_121), xi_124), xi_125), xi_126);
            const __m256 xi_128 = _mm256_mul_ps(xi_127, _mm256_set_ps(omega_even, omega_even, omega_even, omega_even, omega_even, omega_even, omega_even, omega_even));
            const __m256 xi_130 = _mm256_mul_ps(xi_128, _mm256_set_ps(0.01984126984126984f, 0.01984126984126984f, 0.01984126984126984f, 0.01984126984126984f, 0.01984126984126984f, 0.01984126984126984f, 0.01984126984126984f, 0.01984126984126984f));
            const __m256 xi_132 = _mm256_add_ps(xi_130, xi_131);
            const __m256 xi_141 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_113, xi_139), xi_140), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1]));
            const __m256 xi_142 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_117, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_137, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_138, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_141, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_15, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_2, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_4, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0]))), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1])), _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear));
            const __m256 xi_143 = _mm256_mul_ps(xi_142, _mm256_set_ps(0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f, 0.125f));
            const __m256 xi_144 = _mm256_add_ps(_mm256_mul_ps(xi_136, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_143, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)));
            const __m256 xi_168 = _mm256_mul_ps(xi_116, _mm256_set_ps(0.050000000000000003f, 0.050000000000000003f, 0.050000000000000003f, 0.050000000000000003f, 0.050000000000000003f, 0.050000000000000003f, 0.050000000000000003f, 0.050000000000000003f));
            const __m256 xi_170 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0])), _mm256_mul_ps(_mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_ps(_mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_ps(_mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_ps(xi_110, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_118, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_129, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_141, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_147, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(_mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1]))), _mm256_mul_ps(_mm256_set_ps(-2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3 + ctr_0 + 1]))), _mm256_mul_ps(_mm256_set_ps(-2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f, -2.0f), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3 + ctr_0 - 1]))), _mm256_mul_ps(_mm256_mul_ps(rho, xi_111), _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f)));
            const __m256 xi_171 = _mm256_mul_ps(xi_170, _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear));
            const __m256 xi_172 = _mm256_mul_ps(xi_171, _mm256_set_ps(0.041666666666666664f, 0.041666666666666664f, 0.041666666666666664f, 0.041666666666666664f, 0.041666666666666664f, 0.041666666666666664f, 0.041666666666666664f, 0.041666666666666664f));
            const __m256 xi_177 = _mm256_add_ps(xi_172, xi_176);
            const __m256 xi_178 = _mm256_add_ps(_mm256_add_ps(xi_168, xi_169), xi_177);
            const __m256 xi_179 = _mm256_add_ps(_mm256_mul_ps(xi_130, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_131, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)));
            const __m256 xi_180 = _mm256_add_ps(xi_136, xi_143);
            const __m256 xi_181 = _mm256_add_ps(_mm256_mul_ps(xi_172, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_176, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)));
            const __m256 xi_182 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_168, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_169, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_181);
            const __m256 xi_183 = _mm256_mul_ps(xi_123, _mm256_set_ps(0.035714285714285712f, 0.035714285714285712f, 0.035714285714285712f, 0.035714285714285712f, 0.035714285714285712f, 0.035714285714285712f, 0.035714285714285712f, 0.035714285714285712f));
            const __m256 xi_196 = _mm256_mul_ps(xi_123, _mm256_set_ps(0.021428571428571429f, 0.021428571428571429f, 0.021428571428571429f, 0.021428571428571429f, 0.021428571428571429f, 0.021428571428571429f, 0.021428571428571429f, 0.021428571428571429f));
            const __m256 xi_208 = _mm256_mul_ps(xi_142, _mm256_set_ps(0.0625f, 0.0625f, 0.0625f, 0.0625f, 0.0625f, 0.0625f, 0.0625f, 0.0625f));
            const __m256 xi_209 = _mm256_mul_ps(xi_128, _mm256_set_ps(0.013888888888888888f, 0.013888888888888888f, 0.013888888888888888f, 0.013888888888888888f, 0.013888888888888888f, 0.013888888888888888f, 0.013888888888888888f, 0.013888888888888888f));
            const __m256 xi_211 = _mm256_add_ps(_mm256_mul_ps(xi_120, _mm256_set_ps(0.041666666666666664f, 0.041666666666666664f, 0.041666666666666664f, 0.041666666666666664f, 0.041666666666666664f, 0.041666666666666664f, 0.041666666666666664f, 0.041666666666666664f)), _mm256_mul_ps(xi_106, _mm256_set_ps(0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f)));
            const __m256 xi_212 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_171, _mm256_set_ps(0.020833333333333332f, 0.020833333333333332f, 0.020833333333333332f, 0.020833333333333332f, 0.020833333333333332f, 0.020833333333333332f, 0.020833333333333332f, 0.020833333333333332f)), _mm256_mul_ps(xi_175, _mm256_set_ps(0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f))), xi_211);
            const __m256 xi_213 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_166, xi_207), xi_208), xi_209), xi_210), xi_212);
            const __m256 xi_221 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_161, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_165, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_207), xi_208), xi_209), xi_210), xi_212);
            const __m256 xi_225 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(u_2, xi_216), xi_146), xi_155), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0])), _mm256_set_ps(xi_217, xi_217, xi_217, xi_217, xi_217, xi_217, xi_217, xi_217));
            const __m256 xi_232 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_224, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_229, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_222), xi_225), xi_227), xi_231);
            const __m256 xi_238 = _mm256_mul_ps(xi_123, _mm256_set_ps(0.0071428571428571426f, 0.0071428571428571426f, 0.0071428571428571426f, 0.0071428571428571426f, 0.0071428571428571426f, 0.0071428571428571426f, 0.0071428571428571426f, 0.0071428571428571426f));
            const __m256 xi_240 = _mm256_mul_ps(xi_128, _mm256_set_ps(0.003968253968253968f, 0.003968253968253968f, 0.003968253968253968f, 0.003968253968253968f, 0.003968253968253968f, 0.003968253968253968f, 0.003968253968253968f, 0.003968253968253968f));
            const __m256 xi_242 = _mm256_add_ps(_mm256_mul_ps(xi_240, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_241, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)));
            const __m256 xi_243 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_116, _mm256_set_ps(0.025000000000000001f, 0.025000000000000001f, 0.025000000000000001f, 0.025000000000000001f, 0.025000000000000001f, 0.025000000000000001f, 0.025000000000000001f, 0.025000000000000001f)), _mm256_mul_ps(xi_102, _mm256_set_ps(0.050000000000000003f, 0.050000000000000003f, 0.050000000000000003f, 0.050000000000000003f, 0.050000000000000003f, 0.050000000000000003f, 0.050000000000000003f, 0.050000000000000003f))), _mm256_mul_ps(xi_238, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_239, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_181), xi_211), xi_242);
            const __m256 xi_245 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(xi_115, _mm256_set_ps(-0.025000000000000001f, -0.025000000000000001f, -0.025000000000000001f, -0.025000000000000001f, -0.025000000000000001f, -0.025000000000000001f, -0.025000000000000001f, -0.025000000000000001f)), _mm256_set_ps(omega_even, omega_even, omega_even, omega_even, omega_even, omega_even, omega_even, omega_even)), _mm256_mul_ps(_mm256_mul_ps(xi_119, _mm256_set_ps(-0.041666666666666664f, -0.041666666666666664f, -0.041666666666666664f, -0.041666666666666664f, -0.041666666666666664f, -0.041666666666666664f, -0.041666666666666664f, -0.041666666666666664f)), _mm256_set_ps(omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk, omega_bulk))), _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(xi_101, xi_96), _mm256_set_ps(-0.050000000000000003f, -0.050000000000000003f, -0.050000000000000003f, -0.050000000000000003f, -0.050000000000000003f, -0.050000000000000003f, -0.050000000000000003f, -0.050000000000000003f)), _mm256_set_ps(xi_100, xi_100, xi_100, xi_100, xi_100, xi_100, xi_100, xi_100))), _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(xi_103, xi_105), _mm256_set_ps(-0.083333333333333329f, -0.083333333333333329f, -0.083333333333333329f, -0.083333333333333329f, -0.083333333333333329f, -0.083333333333333329f, -0.083333333333333329f, -0.083333333333333329f)), _mm256_set_ps(xi_104, xi_104, xi_104, xi_104, xi_104, xi_104, xi_104, xi_104))), xi_177), xi_238), xi_239), xi_240), xi_241);
            const __m256 xi_247 = _mm256_mul_ps(xi_208, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f));
            const __m256 xi_248 = _mm256_mul_ps(xi_123, _mm256_set_ps(0.017857142857142856f, 0.017857142857142856f, 0.017857142857142856f, 0.017857142857142856f, 0.017857142857142856f, 0.017857142857142856f, 0.017857142857142856f, 0.017857142857142856f));
            const __m256 xi_250 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(xi_205, xi_212), xi_242), xi_246), xi_247), xi_248), xi_249);
            const __m256 xi_252 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(rho, u_0), u_2), xi_137), xi_185), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1])), _mm256_set_ps(xi_217, xi_217, xi_217, xi_217, xi_217, xi_217, xi_217, xi_217));
            const __m256 xi_253 = _mm256_add_ps(_mm256_mul_ps(xi_251, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_252, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)));
            const __m256 xi_259 = _mm256_add_ps(xi_251, xi_252);
            const __m256 xi_261 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_227, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_231, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_222), xi_224), xi_225), xi_229);
            const __m256 xi_262 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_202, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_204, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_212), xi_242), xi_246), xi_247), xi_248), xi_249);
            const __m256 forceTerm_0 = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_17, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_18, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_19, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_17, _mm256_set_ps(xi_20, xi_20, xi_20, xi_20, xi_20, xi_20, xi_20, xi_20))), _mm256_mul_ps(xi_18, _mm256_set_ps(xi_20, xi_20, xi_20, xi_20, xi_20, xi_20, xi_20, xi_20))), _mm256_mul_ps(xi_19, _mm256_set_ps(xi_20, xi_20, xi_20, xi_20, xi_20, xi_20, xi_20, xi_20)));
            const __m256 forceTerm_1 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_23, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_21), xi_32);
            const __m256 forceTerm_2 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_21, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_23), xi_32);
            const __m256 forceTerm_3 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_33, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_38, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_35);
            const __m256 forceTerm_4 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_35, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_38, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), xi_33);
            const __m256 forceTerm_5 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_41, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_39), xi_44);
            const __m256 forceTerm_6 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_39, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), xi_41), xi_44);
            const __m256 forceTerm_7 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_49, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_57, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_67, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)));
            const __m256 forceTerm_8 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_57, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_68, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_69, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)));
            const __m256 forceTerm_9 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_49, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_69, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_71, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)));
            const __m256 forceTerm_10 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_67, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_68, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_71, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)));
            const __m256 forceTerm_11 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_74, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_76, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_81, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)));
            const __m256 forceTerm_12 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_74, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_82, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_83, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)));
            const __m256 forceTerm_13 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_49, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_84, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_91, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)));
            const __m256 forceTerm_14 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_68, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_84, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_92, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)));
            const __m256 forceTerm_15 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_76, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_83, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_93, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)));
            const __m256 forceTerm_16 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_81, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_82, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_93, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)));
            const __m256 forceTerm_17 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_49, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_92, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_94, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)));
            const __m256 forceTerm_18 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_68, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_91, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_94, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)));
            _mm256_store_ps(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_128, _mm256_set_ps(0.023809523809523808f, 0.023809523809523808f, 0.023809523809523808f, 0.023809523809523808f, 0.023809523809523808f, 0.023809523809523808f, 0.023809523809523808f, 0.023809523809523808f)), _mm256_mul_ps(xi_99, _mm256_set_ps(0.14285714285714285f, 0.14285714285714285f, 0.14285714285714285f, 0.14285714285714285f, 0.14285714285714285f, 0.14285714285714285f, 0.14285714285714285f, 0.14285714285714285f))), _mm256_mul_ps(xi_123, _mm256_set_ps(0.042857142857142858f, 0.042857142857142858f, 0.042857142857142858f, 0.042857142857142858f, 0.042857142857142858f, 0.042857142857142858f, 0.042857142857142858f, 0.042857142857142858f))), _mm256_mul_ps(xi_109, _mm256_set_ps(0.085714285714285715f, 0.085714285714285715f, 0.085714285714285715f, 0.085714285714285715f, 0.085714285714285715f, 0.085714285714285715f, 0.085714285714285715f, 0.085714285714285715f))), _mm256_mul_ps(xi_116, _mm256_set_ps(0.10000000000000001f, 0.10000000000000001f, 0.10000000000000001f, 0.10000000000000001f, 0.10000000000000001f, 0.10000000000000001f, 0.10000000000000001f, 0.10000000000000001f))), _mm256_mul_ps(xi_102, _mm256_set_ps(0.20000000000000001f, 0.20000000000000001f, 0.20000000000000001f, 0.20000000000000001f, 0.20000000000000001f, 0.20000000000000001f, 0.20000000000000001f, 0.20000000000000001f))), _mm256_mul_ps(xi_120, _mm256_set_ps(-0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f))), _mm256_mul_ps(xi_106, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), forceTerm_0), _mm256_load_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + ctr_0])));
            _mm256_storeu_ps(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + _stride_pdfs_tmp_3 + ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_129, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_132, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_144, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_167, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_178, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(_mm256_mul_ps(xi_122, _mm256_set_ps(0.014285714285714285f, 0.014285714285714285f, 0.014285714285714285f, 0.014285714285714285f, 0.014285714285714285f, 0.014285714285714285f, 0.014285714285714285f, 0.014285714285714285f)), _mm256_set_ps(omega_even, omega_even, omega_even, omega_even, omega_even, omega_even, omega_even, omega_even))), _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(xi_108, xi_96), _mm256_set_ps(0.028571428571428571f, 0.028571428571428571f, 0.028571428571428571f, 0.028571428571428571f, 0.028571428571428571f, 0.028571428571428571f, 0.028571428571428571f, 0.028571428571428571f)), _mm256_set_ps(xi_107, xi_107, xi_107, xi_107, xi_107, xi_107, xi_107, xi_107))), forceTerm_1));
            _mm256_storeu_ps(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 2 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_123, _mm256_set_ps(0.014285714285714285f, 0.014285714285714285f, 0.014285714285714285f, 0.014285714285714285f, 0.014285714285714285f, 0.014285714285714285f, 0.014285714285714285f, 0.014285714285714285f)), _mm256_mul_ps(xi_109, _mm256_set_ps(0.028571428571428571f, 0.028571428571428571f, 0.028571428571428571f, 0.028571428571428571f, 0.028571428571428571f, 0.028571428571428571f, 0.028571428571428571f, 0.028571428571428571f))), forceTerm_2), xi_167), xi_179), xi_180), xi_182), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0])));
            _mm256_storeu_ps(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 3 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_171, _mm256_set_ps(0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f)), _mm256_mul_ps(xi_175, _mm256_set_ps(0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f))), _mm256_mul_ps(xi_183, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_184, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), forceTerm_3), xi_179), xi_195), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3 + ctr_0 + 1])));
            _mm256_storeu_ps(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 4 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_132, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_183, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_184, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_189, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_195, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(_mm256_mul_ps(xi_170, _mm256_set_ps(0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f, 0.083333333333333329f)), _mm256_set_ps(omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear, omega_shear))), _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(xi_134, xi_174), _mm256_set_ps(0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f, 0.33333333333333331f)), _mm256_set_ps(xi_173, xi_173, xi_173, xi_173, xi_173, xi_173, xi_173, xi_173))), forceTerm_4));
            _mm256_storeu_ps(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 5 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_139, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_178, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_180, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_196, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_197, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_206, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(_mm256_mul_ps(xi_127, _mm256_set_ps(0.015873015873015872f, 0.015873015873015872f, 0.015873015873015872f, 0.015873015873015872f, 0.015873015873015872f, 0.015873015873015872f, 0.015873015873015872f, 0.015873015873015872f)), _mm256_set_ps(omega_even, omega_even, omega_even, omega_even, omega_even, omega_even, omega_even, omega_even))), _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(xi_96, xi_98), _mm256_set_ps(0.095238095238095233f, 0.095238095238095233f, 0.095238095238095233f, 0.095238095238095233f, 0.095238095238095233f, 0.095238095238095233f, 0.095238095238095233f, 0.095238095238095233f)), _mm256_set_ps(xi_97, xi_97, xi_97, xi_97, xi_97, xi_97, xi_97, xi_97))), forceTerm_5));
            _mm256_storeu_ps(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 6 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_128, _mm256_set_ps(0.015873015873015872f, 0.015873015873015872f, 0.015873015873015872f, 0.015873015873015872f, 0.015873015873015872f, 0.015873015873015872f, 0.015873015873015872f, 0.015873015873015872f)), _mm256_mul_ps(xi_99, _mm256_set_ps(0.095238095238095233f, 0.095238095238095233f, 0.095238095238095233f, 0.095238095238095233f, 0.095238095238095233f, 0.095238095238095233f, 0.095238095238095233f, 0.095238095238095233f))), _mm256_mul_ps(xi_196, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_197, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), forceTerm_6), xi_144), xi_182), xi_206), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0])));
            _mm256_storeu_ps(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 7 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_7, xi_213), xi_214), xi_219), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3 + ctr_0 + 1])));
            _mm256_store_ps(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 8 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_8, xi_194), xi_213), xi_220), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3 + ctr_0 - 1])));
            _mm256_storeu_ps(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 9 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_9, xi_214), xi_220), xi_221), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1])));
            _mm256_storeu_ps(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 10 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_10, xi_194), xi_219), xi_221), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3 + ctr_0 - 1])));
            _mm256_storeu_ps(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 11 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_11, xi_232), xi_237), xi_243), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0])));
            _mm256_storeu_ps(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 12 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_156, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_232, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_244, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_245, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), forceTerm_12));
            _mm256_storeu_ps(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 13 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_13, xi_250), xi_253), xi_258), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1])));
            _mm256_storeu_ps(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 14 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_14, xi_250), xi_259), xi_260), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3 + ctr_0 - 1])));
            _mm256_storeu_ps(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 15 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(xi_145, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f)), _mm256_mul_ps(xi_237, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_245, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), _mm256_mul_ps(xi_261, _mm256_set_ps(-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f))), forceTerm_15));
            _mm256_store_ps(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 16 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_16, xi_243), xi_244), xi_261), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0])));
            _mm256_storeu_ps(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 17 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_17, xi_258), xi_259), xi_262), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1])));
            _mm256_storeu_ps(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 18 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(forceTerm_18, xi_253), xi_260), xi_262), _mm256_loadu_ps(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1])));
          }
          for (int64_t ctr_0 = (int64_t)((_size_force_0 - 2) / (8)) * (8) + 1; ctr_0 < _size_force_0 - 1; ctr_0 += 1) {

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
            const float xi_2 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3 + ctr_0 - 1] + _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3 + ctr_0 - 1];
            const float xi_3 = xi_2 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1] + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3 + ctr_0 - 1];
            const float xi_4 = _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3 + ctr_0 + 1] + _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3 + ctr_0];
            const float xi_5 = _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0] + _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0];
            const float xi_6 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0] + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3 + ctr_0];
            const float xi_7 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3 + ctr_0 + 1] + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1];
            const float xi_8 = xi_7 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1];
            const float xi_9 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0] + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0];
            const float xi_11 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1] + _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3 + ctr_0 + 1];
            const float xi_12 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3 + ctr_0 - 1] - _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3 + ctr_0 - 1];
            const float xi_13 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1] + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0];
            const float xi_14 = -_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0];
            const float xi_15 = -_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3 + ctr_0 - 1];
            const float xi_16 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0] + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1];
            const float xi_21 = 0.16666666666666666f * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0];
            const float xi_22 = 0.083333333333333329f * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0];
            const float xi_33 = 0.16666666666666666f * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0];
            const float xi_34 = 0.083333333333333329f * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0];
            const float xi_39 = 0.16666666666666666f * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0];
            const float xi_40 = 0.083333333333333329f * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0];
            const float xi_58 = 0.25f * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0];
            const float xi_63 = xi_62 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0];
            const float xi_98 = random_3_0 - 0.5f;
            const float xi_101 = random_3_2 - 0.5f;
            const float xi_103 = random_1_1 - 0.5f;
            const float xi_108 = random_3_1 - 0.5f;
            const float xi_112 = -_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + ctr_0];
            const float xi_118 = xi_11 + xi_3;
            const float xi_121 = 2.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0] + 2.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0] + 2.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0] + 2.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0] + 5.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3 + ctr_0 + 1] + 5.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3 + ctr_0 - 1];
            const float xi_124 = 2.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1];
            const float xi_125 = 2.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3 + ctr_0 - 1];
            const float xi_126 = 2.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1] + 2.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1];
            const float xi_129 = -_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3 + ctr_0];
            const float xi_133 = random_0_1 - 0.5f;
            const float xi_137 = -_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1];
            const float xi_138 = -_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1];
            const float xi_139 = -_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3 + ctr_0];
            const float xi_140 = -_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0];
            const float xi_145 = -_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0];
            const float xi_146 = xi_14 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0];
            const float xi_147 = -_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0];
            const float xi_148 = xi_147 + _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3 + ctr_0];
            const float xi_149 = xi_145 + xi_146 + xi_148 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0];
            const float xi_151 = random_2_3 - 0.5f;
            const float xi_155 = -_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0];
            const float xi_156 = -_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0];
            const float xi_157 = 2.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3 + ctr_0 - 1];
            const float xi_158 = -2.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3 + ctr_0 + 1] + 2.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3 + ctr_0 - 1];
            const float xi_159 = xi_148 + xi_155 + xi_156 - xi_157 + xi_158 + xi_5 + 2.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1];
            const float xi_161 = xi_159 * xi_160;
            const float xi_162 = random_1_2 - 0.5f;
            const float xi_174 = random_0_0 - 0.5f;
            const float xi_185 = xi_15 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1];
            const float xi_186 = xi_138 + xi_185;
            const float xi_187 = xi_186 - _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3 + ctr_0 + 1] + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3 + ctr_0 - 1] + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1];
            const float xi_188 = random_2_1 - 0.5f;
            const float xi_189 = -_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3 + ctr_0 - 1];
            const float xi_190 = -xi_157 - xi_158 - xi_186 - xi_189 - xi_7 + 2.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1];
            const float xi_191 = xi_160 * xi_190;
            const float xi_192 = random_2_0 - 0.5f;
            const float xi_198 = xi_145 + xi_155 + _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0];
            const float xi_199 = -xi_139 - xi_198 - _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0] - _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0];
            const float xi_200 = random_2_2 - 0.5f;
            const float xi_201 = -xi_124 - xi_125 + xi_126 + xi_140 + xi_198 + xi_6;
            const float xi_202 = xi_160 * xi_201;
            const float xi_203 = random_1_3 - 0.5f;
            const float xi_224 = xi_159 * xi_223;
            const float xi_227 = xi_149 * xi_226;
            const float xi_233 = xi_201 * xi_223;
            const float xi_234 = xi_199 * xi_226;
            const float xi_254 = xi_187 * xi_226;
            const float xi_255 = xi_190 * xi_223;
            const float xi_23 = rr_0 * xi_22;
            const float xi_35 = rr_0 * xi_34;
            const float xi_41 = rr_0 * xi_40;
            const float xi_46 = xi_45 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0];
            const float xi_50 = xi_45 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0];
            const float xi_72 = xi_45 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0];
            const float vel0Term = xi_3 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3 + ctr_0 - 1];
            const float vel1Term = xi_4 + xi_5;
            const float vel2Term = xi_6 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1];
            const float delta_rho = vel0Term + vel1Term + vel2Term + xi_8 + xi_9 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0] + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + ctr_0];
            const float rho = delta_rho + 1.0f;
            const float xi_95 = kT * rho;
            const float xi_96 = powf(xi_95 * (1.0f - ((-omega_even + 1.0f) * (-omega_even + 1.0f))), 0.5f);
            const float xi_99 = xi_96 * xi_97 * xi_98;
            const float xi_102 = xi_100 * xi_101 * xi_96;
            const float xi_105 = powf(xi_95 * (1.0f - ((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f))), 0.5f);
            const float xi_106 = xi_103 * xi_104 * xi_105;
            const float xi_109 = xi_107 * xi_108 * xi_96;
            const float xi_131 = xi_99 * 0.11904761904761904f;
            const float xi_134 = powf(xi_95 * (1.0f - ((-omega_shear + 1.0f) * (-omega_shear + 1.0f))), 0.5f);
            const float xi_135 = xi_134 * 0.5f;
            const float xi_136 = xi_133 * xi_135;
            const float xi_152 = powf(xi_95 * (1.0f - ((-omega_odd + 1.0f) * (-omega_odd + 1.0f))), 0.5f);
            const float xi_153 = xi_152 * 1.4142135623730951f;
            const float xi_154 = xi_153 * 0.5f;
            const float xi_163 = xi_104 * xi_152;
            const float xi_164 = xi_163 * 0.16666666666666666f;
            const float xi_165 = xi_162 * xi_164;
            const float xi_166 = xi_161 + xi_165;
            const float xi_167 = xi_149 * xi_150 + xi_151 * xi_154 + xi_166;
            const float xi_169 = xi_102 * 0.10000000000000001f;
            const float xi_175 = xi_134 * xi_173 * xi_174;
            const float xi_176 = xi_175 * 0.16666666666666666f;
            const float xi_184 = xi_109 * 0.071428571428571425f;
            const float xi_193 = xi_164 * xi_192;
            const float xi_194 = xi_191 + xi_193;
            const float xi_195 = xi_150 * xi_187 + xi_154 * xi_188 + xi_194;
            const float xi_197 = xi_109 * 0.042857142857142858f;
            const float xi_204 = xi_164 * xi_203;
            const float xi_205 = xi_202 + xi_204;
            const float xi_206 = xi_150 * xi_199 + xi_154 * xi_200 + xi_205;
            const float xi_207 = xi_133 * xi_134 * 0.25f;
            const float xi_210 = xi_99 * 0.083333333333333329f;
            const float xi_214 = -xi_191 - xi_193;
            const float xi_215 = xi_135 * (random_0_2 - 0.5f);
            const float xi_222 = xi_135 * (random_1_0 - 0.5f);
            const float xi_228 = xi_163 * 0.083333333333333329f;
            const float xi_229 = xi_162 * xi_228;
            const float xi_230 = xi_153 * 0.25f;
            const float xi_231 = xi_151 * xi_230;
            const float xi_235 = xi_203 * xi_228;
            const float xi_236 = xi_200 * xi_230;
            const float xi_237 = -xi_233 + xi_234 - xi_235 + xi_236;
            const float xi_239 = xi_109 * 0.014285714285714285f;
            const float xi_241 = xi_99 * 0.023809523809523808f;
            const float xi_244 = xi_233 - xi_234 + xi_235 - xi_236;
            const float xi_246 = -xi_207;
            const float xi_249 = xi_109 * 0.035714285714285712f;
            const float xi_251 = xi_135 * (random_0_3 - 0.5f);
            const float xi_256 = xi_188 * xi_230;
            const float xi_257 = xi_192 * xi_228;
            const float xi_258 = -xi_254 + xi_255 - xi_256 + xi_257;
            const float xi_260 = xi_254 - xi_255 + xi_256 - xi_257;
            const float xi_0 = ((1.0f) / (rho));
            const float xi_10 = xi_0 * 0.5f;
            const float u_0 = xi_0 * (vel0Term - xi_11 - xi_8) + xi_10 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0];
            const float xi_17 = u_0 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0];
            const float xi_28 = xi_17 * 0.16666666666666666f;
            const float xi_29 = -xi_28;
            const float xi_30 = xi_17 * 0.083333333333333329f;
            const float xi_31 = omega_shear * xi_30 + xi_29;
            const float xi_48 = xi_17 * xi_47 + xi_29;
            const float xi_49 = xi_34 - xi_46 + xi_48;
            const float xi_52 = xi_17 * xi_51;
            const float xi_59 = u_0 * xi_58;
            const float xi_64 = u_0 * xi_63;
            const float xi_68 = -xi_34 + xi_46 + xi_48;
            const float xi_75 = omega_shear * u_0 * -0.083333333333333329f * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0];
            const float xi_85 = u_0 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0];
            const float xi_86 = xi_85 * 0.25f;
            const float xi_89 = xi_62 * xi_85;
            const float xi_111 = (u_0 * u_0);
            const float u_1 = xi_0 * (vel1Term - xi_12 - xi_13 - xi_9) + xi_10 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0];
            const float xi_18 = u_1 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0];
            const float xi_26 = xi_18 * 0.16666666666666666f;
            const float xi_36 = omega_shear * u_1 * -0.083333333333333329f * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0];
            const float xi_42 = -xi_26;
            const float xi_43 = xi_18 * 0.083333333333333329f;
            const float xi_53 = xi_18 * xi_51;
            const float xi_60 = u_1 * 0.25f;
            const float xi_61 = xi_60 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0];
            const float xi_65 = u_1 * xi_62;
            const float xi_66 = xi_65 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0];
            const float xi_67 = xi_59 + xi_61 - xi_64 - xi_66;
            const float xi_69 = -xi_59 - xi_61 + xi_64 + xi_66;
            const float xi_77 = xi_60 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0];
            const float xi_79 = xi_65 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0];
            const float xi_110 = rho * (u_1 * u_1);
            const float xi_117 = -xi_110;
            const float xi_216 = rho * u_1;
            const float xi_218 = xi_217 * (u_0 * xi_216 + xi_12 - _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1] + _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3 + ctr_0 + 1]);
            const float xi_219 = -xi_215 - xi_218;
            const float xi_220 = xi_215 + xi_218;
            const float u_2 = xi_0 * (vel2Term - xi_14 - xi_15 - xi_16 - _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1] - _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0] - _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0]) + xi_10 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0];
            const float xi_19 = u_2 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0];
            const float xi_24 = xi_19 * 0.16666666666666666f;
            const float xi_25 = -xi_24;
            const float xi_27 = xi_19 * 0.083333333333333329f;
            const float xi_32 = -omega_shear * xi_26 + omega_shear * xi_27 + xi_18 * 0.33333333333333331f + xi_25 + xi_31;
            const float xi_37 = omega_shear * u_2 * -0.083333333333333329f * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0];
            const float xi_38 = omega_shear * xi_28 + u_0 * -0.33333333333333331f * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0] + xi_24 + xi_26 + xi_36 + xi_37;
            const float xi_44 = -omega_shear * xi_24 + omega_shear * xi_43 + xi_19 * 0.33333333333333331f + xi_31 + xi_42;
            const float xi_54 = xi_19 * xi_51;
            const float xi_55 = xi_18 * xi_47 + xi_42 + xi_52 + xi_53 + xi_54;
            const float xi_56 = -xi_22 + xi_50 + xi_55;
            const float xi_57 = xi_27 + xi_37 + xi_56;
            const float xi_70 = xi_22 - xi_50 + xi_55;
            const float xi_71 = xi_27 + xi_37 + xi_70;
            const float xi_73 = xi_19 * xi_47 + xi_25;
            const float xi_74 = -xi_40 + xi_72 + xi_73;
            const float xi_76 = xi_30 + xi_56 + xi_75;
            const float xi_78 = u_2 * xi_58;
            const float xi_80 = u_2 * xi_63;
            const float xi_81 = -xi_77 - xi_78 + xi_79 + xi_80;
            const float xi_82 = xi_30 + xi_70 + xi_75;
            const float xi_83 = xi_77 + xi_78 - xi_79 - xi_80;
            const float xi_84 = xi_36 + xi_43 + xi_52 + xi_53 + xi_54 + xi_74;
            const float xi_87 = u_2 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0];
            const float xi_88 = xi_87 * 0.25f;
            const float xi_90 = xi_62 * xi_87;
            const float xi_91 = xi_86 + xi_88 - xi_89 - xi_90;
            const float xi_92 = -xi_86 - xi_88 + xi_89 + xi_90;
            const float xi_93 = xi_40 - xi_72 + xi_73;
            const float xi_94 = xi_36 + xi_43 + xi_52 + xi_53 + xi_54 + xi_93;
            const float xi_113 = rho * (u_2 * u_2);
            const float xi_114 = xi_112 + xi_113 * 0.66666666666666663f + 3.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0] + 3.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3 + ctr_0];
            const float xi_115 = rho * xi_111 * 1.6666666666666667f + xi_110 * 0.66666666666666663f + xi_114 - 3.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0] - 3.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0] - 3.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0] - 3.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0] + 3.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0] + 3.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3 + ctr_0];
            const float xi_116 = omega_even * xi_115;
            const float xi_119 = rho * xi_111 - xi_112 + xi_113 - xi_117 - xi_118 - xi_13 - xi_16 - xi_5;
            const float xi_120 = omega_bulk * xi_119;
            const float xi_122 = xi_110 * 2.3333333333333335f + xi_114 + xi_121 - 2.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0] - 2.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3 + ctr_0] - 5.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1] - 5.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1] - 5.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1] - 5.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3 + ctr_0 - 1];
            const float xi_123 = omega_even * xi_122;
            const float xi_127 = xi_112 + xi_113 * 3.0f + xi_121 + xi_124 + xi_125 + xi_126 - 4.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0] - 4.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3 + ctr_0] - 7.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3 + ctr_0 - 1] - 7.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1] - 7.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3 + ctr_0 + 1] - 7.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3 + ctr_0 - 1] + 5.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0] + 5.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3 + ctr_0];
            const float xi_128 = omega_even * xi_127;
            const float xi_130 = xi_128 * 0.01984126984126984f;
            const float xi_132 = xi_130 + xi_131;
            const float xi_141 = xi_113 + xi_139 + xi_140 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1];
            const float xi_142 = omega_shear * (-xi_117 - xi_137 - xi_138 - xi_141 - xi_15 - xi_2 - xi_4 - _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0] + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1]);
            const float xi_143 = xi_142 * 0.125f;
            const float xi_144 = -xi_136 - xi_143;
            const float xi_168 = xi_116 * 0.050000000000000003f;
            const float xi_170 = rho * xi_111 * 2.0f - xi_110 - xi_118 - xi_129 - xi_141 - xi_147 - 2.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3 + ctr_0 + 1] - 2.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3 + ctr_0 - 1] + 2.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0] + 2.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0] + 2.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0] + 2.0f * _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0] - _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1];
            const float xi_171 = omega_shear * xi_170;
            const float xi_172 = xi_171 * 0.041666666666666664f;
            const float xi_177 = xi_172 + xi_176;
            const float xi_178 = xi_168 + xi_169 + xi_177;
            const float xi_179 = -xi_130 - xi_131;
            const float xi_180 = xi_136 + xi_143;
            const float xi_181 = -xi_172 - xi_176;
            const float xi_182 = -xi_168 - xi_169 + xi_181;
            const float xi_183 = xi_123 * 0.035714285714285712f;
            const float xi_196 = xi_123 * 0.021428571428571429f;
            const float xi_208 = xi_142 * 0.0625f;
            const float xi_209 = xi_128 * 0.013888888888888888f;
            const float xi_211 = xi_106 * 0.083333333333333329f + xi_120 * 0.041666666666666664f;
            const float xi_212 = xi_171 * 0.020833333333333332f + xi_175 * 0.083333333333333329f + xi_211;
            const float xi_213 = xi_166 + xi_207 + xi_208 + xi_209 + xi_210 + xi_212;
            const float xi_221 = -xi_161 - xi_165 + xi_207 + xi_208 + xi_209 + xi_210 + xi_212;
            const float xi_225 = xi_217 * (u_2 * xi_216 + xi_146 + xi_155 + _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0]);
            const float xi_232 = xi_222 - xi_224 + xi_225 + xi_227 - xi_229 + xi_231;
            const float xi_238 = xi_123 * 0.0071428571428571426f;
            const float xi_240 = xi_128 * 0.003968253968253968f;
            const float xi_242 = -xi_240 - xi_241;
            const float xi_243 = xi_102 * 0.050000000000000003f + xi_116 * 0.025000000000000001f + xi_181 + xi_211 - xi_238 - xi_239 + xi_242;
            const float xi_245 = omega_bulk * xi_119 * -0.041666666666666664f + omega_even * xi_115 * -0.025000000000000001f + xi_100 * xi_101 * xi_96 * -0.050000000000000003f + xi_103 * xi_104 * xi_105 * -0.083333333333333329f + xi_177 + xi_238 + xi_239 + xi_240 + xi_241;
            const float xi_247 = -xi_208;
            const float xi_248 = xi_123 * 0.017857142857142856f;
            const float xi_250 = xi_205 + xi_212 + xi_242 + xi_246 + xi_247 + xi_248 + xi_249;
            const float xi_252 = xi_217 * (rho * u_0 * u_2 + xi_137 + xi_185 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1]);
            const float xi_253 = -xi_251 - xi_252;
            const float xi_259 = xi_251 + xi_252;
            const float xi_261 = xi_222 + xi_224 + xi_225 - xi_227 + xi_229 - xi_231;
            const float xi_262 = -xi_202 - xi_204 + xi_212 + xi_242 + xi_246 + xi_247 + xi_248 + xi_249;
            const float forceTerm_0 = xi_17 * xi_20 - xi_17 + xi_18 * xi_20 - xi_18 + xi_19 * xi_20 - xi_19;
            const float forceTerm_1 = xi_21 - xi_23 + xi_32;
            const float forceTerm_2 = -xi_21 + xi_23 + xi_32;
            const float forceTerm_3 = -xi_33 + xi_35 - xi_38;
            const float forceTerm_4 = xi_33 - xi_35 - xi_38;
            const float forceTerm_5 = xi_39 - xi_41 + xi_44;
            const float forceTerm_6 = -xi_39 + xi_41 + xi_44;
            const float forceTerm_7 = -xi_49 - xi_57 - xi_67;
            const float forceTerm_8 = -xi_57 - xi_68 - xi_69;
            const float forceTerm_9 = -xi_49 - xi_69 - xi_71;
            const float forceTerm_10 = -xi_67 - xi_68 - xi_71;
            const float forceTerm_11 = -xi_74 - xi_76 - xi_81;
            const float forceTerm_12 = -xi_74 - xi_82 - xi_83;
            const float forceTerm_13 = -xi_49 - xi_84 - xi_91;
            const float forceTerm_14 = -xi_68 - xi_84 - xi_92;
            const float forceTerm_15 = -xi_76 - xi_83 - xi_93;
            const float forceTerm_16 = -xi_81 - xi_82 - xi_93;
            const float forceTerm_17 = -xi_49 - xi_92 - xi_94;
            const float forceTerm_18 = -xi_68 - xi_91 - xi_94;
            _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + ctr_0] = forceTerm_0 + xi_102 * 0.20000000000000001f - xi_106 + xi_109 * 0.085714285714285715f + xi_116 * 0.10000000000000001f + xi_120 * -0.5f + xi_123 * 0.042857142857142858f + xi_128 * 0.023809523809523808f + xi_99 * 0.14285714285714285f + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + ctr_0];
            _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + _stride_pdfs_tmp_3 + ctr_0] = forceTerm_1 + omega_even * xi_122 * 0.014285714285714285f + xi_107 * xi_108 * xi_96 * 0.028571428571428571f - xi_129 - xi_132 - xi_144 - xi_167 - xi_178;
            _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 2 * _stride_pdfs_tmp_3 + ctr_0] = forceTerm_2 + xi_109 * 0.028571428571428571f + xi_123 * 0.014285714285714285f + xi_167 + xi_179 + xi_180 + xi_182 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0];
            _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 3 * _stride_pdfs_tmp_3 + ctr_0] = forceTerm_3 + xi_171 * 0.083333333333333329f + xi_175 * 0.33333333333333331f + xi_179 - xi_183 - xi_184 + xi_195 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3 + ctr_0 + 1];
            _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 4 * _stride_pdfs_tmp_3 + ctr_0] = forceTerm_4 + omega_shear * xi_170 * 0.083333333333333329f - xi_132 + xi_134 * xi_173 * xi_174 * 0.33333333333333331f - xi_183 - xi_184 - xi_189 - xi_195;
            _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 5 * _stride_pdfs_tmp_3 + ctr_0] = forceTerm_5 + omega_even * xi_127 * 0.015873015873015872f - xi_139 - xi_178 - xi_180 - xi_196 - xi_197 - xi_206 + xi_96 * xi_97 * xi_98 * 0.095238095238095233f;
            _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 6 * _stride_pdfs_tmp_3 + ctr_0] = forceTerm_6 + xi_128 * 0.015873015873015872f + xi_144 + xi_182 - xi_196 - xi_197 + xi_206 + xi_99 * 0.095238095238095233f + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0];
            _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 7 * _stride_pdfs_tmp_3 + ctr_0] = forceTerm_7 + xi_213 + xi_214 + xi_219 + _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3 + ctr_0 + 1];
            _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 8 * _stride_pdfs_tmp_3 + ctr_0] = forceTerm_8 + xi_194 + xi_213 + xi_220 + _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3 + ctr_0 - 1];
            _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 9 * _stride_pdfs_tmp_3 + ctr_0] = forceTerm_9 + xi_214 + xi_220 + xi_221 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1];
            _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 10 * _stride_pdfs_tmp_3 + ctr_0] = forceTerm_10 + xi_194 + xi_219 + xi_221 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3 + ctr_0 - 1];
            _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 11 * _stride_pdfs_tmp_3 + ctr_0] = forceTerm_11 + xi_232 + xi_237 + xi_243 + _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0];
            _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 12 * _stride_pdfs_tmp_3 + ctr_0] = forceTerm_12 - xi_156 - xi_232 - xi_244 - xi_245;
            _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 13 * _stride_pdfs_tmp_3 + ctr_0] = forceTerm_13 + xi_250 + xi_253 + xi_258 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1];
            _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 14 * _stride_pdfs_tmp_3 + ctr_0] = forceTerm_14 + xi_250 + xi_259 + xi_260 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3 + ctr_0 - 1];
            _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 15 * _stride_pdfs_tmp_3 + ctr_0] = forceTerm_15 - xi_145 - xi_237 - xi_245 - xi_261;
            _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 16 * _stride_pdfs_tmp_3 + ctr_0] = forceTerm_16 + xi_243 + xi_244 + xi_261 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0];
            _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 17 * _stride_pdfs_tmp_3 + ctr_0] = forceTerm_17 + xi_258 + xi_259 + xi_262 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1];
            _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 18 * _stride_pdfs_tmp_3 + ctr_0] = forceTerm_18 + xi_253 + xi_260 + xi_262 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1];
          }
        }
      }
    }
  }
}
} // namespace internal_482af0085ade5d86c474548ac4b32084

void StreamCollideSweepThermalizedSinglePrecisionAVX::run(IBlock *block) {
  if (!this->configured_)
    WALBERLA_ABORT("This Sweep contains a configure function that needs to be called manually")

  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);
  auto force = block->getData<field::GhostLayerField<float, 3>>(forceID);
  field::GhostLayerField<float, 19> *pdfs_tmp;
  {
    if (cache_pdfs_.find(block) == cache_pdfs_.end()) {
      pdfs_tmp = pdfs->cloneUninitialized();
      cache_pdfs_[block] = pdfs_tmp;
    } else {
      pdfs_tmp = cache_pdfs_[block];
    }
  }

  auto &omega_shear = this->omega_shear_;
  auto &kT = this->kT_;
  auto &block_offset_1 = this->block_offset_1_;
  auto &omega_bulk = this->omega_bulk_;
  auto &time_step = this->time_step_;
  auto &block_offset_2 = this->block_offset_2_;
  auto &block_offset_0 = this->block_offset_0_;
  auto &omega_odd = this->omega_odd_;
  auto &seed = this->seed_;
  auto &omega_even = this->omega_even_;
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force->nrOfGhostLayers()))
  float *RESTRICT const _data_force = force->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs->nrOfGhostLayers()))
  float *RESTRICT const _data_pdfs = pdfs->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)pdfs->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_tmp->nrOfGhostLayers()))
  float *RESTRICT _data_pdfs_tmp = pdfs_tmp->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs_tmp->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)pdfs_tmp->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(), int64_t(int64_c(force->xSize()) + 2))
  const int64_t _size_force_0 = int64_t(int64_c(force->xSize()) + 2);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(), int64_t(int64_c(force->ySize()) + 2))
  const int64_t _size_force_1 = int64_t(int64_c(force->ySize()) + 2);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(), int64_t(int64_c(force->zSize()) + 2))
  const int64_t _size_force_2 = int64_t(int64_c(force->zSize()) + 2);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0)
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  const int64_t _stride_pdfs_tmp_1 = int64_t(pdfs_tmp->yStride());
  const int64_t _stride_pdfs_tmp_2 = int64_t(pdfs_tmp->zStride());
  const int64_t _stride_pdfs_tmp_3 = int64_t(1 * int64_t(pdfs_tmp->fStride()));
  internal_482af0085ade5d86c474548ac4b32084::streamcollidesweepthermalizedsingleprecisionavx_streamcollidesweepthermalizedsingleprecisionavx(_data_force, _data_pdfs, _data_pdfs_tmp, _size_force_0, _size_force_1, _size_force_2, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_pdfs_tmp_3, block_offset_0, block_offset_1, block_offset_2, kT, omega_bulk, omega_even, omega_odd, omega_shear, seed, time_step);
  pdfs->swapDataPointers(pdfs_tmp);
}

void StreamCollideSweepThermalizedSinglePrecisionAVX::runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks, const CellInterval &globalCellInterval, cell_idx_t ghostLayers, IBlock *block) {
  if (!this->configured_)
    WALBERLA_ABORT("This Sweep contains a configure function that needs to be called manually")

  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);
  auto force = block->getData<field::GhostLayerField<float, 3>>(forceID);
  field::GhostLayerField<float, 19> *pdfs_tmp;
  {
    if (cache_pdfs_.find(block) == cache_pdfs_.end()) {
      pdfs_tmp = pdfs->cloneUninitialized();
      cache_pdfs_[block] = pdfs_tmp;
    } else {
      pdfs_tmp = cache_pdfs_[block];
    }
  }

  auto &omega_shear = this->omega_shear_;
  auto &kT = this->kT_;
  auto &block_offset_1 = this->block_offset_1_;
  auto &omega_bulk = this->omega_bulk_;
  auto &time_step = this->time_step_;
  auto &block_offset_2 = this->block_offset_2_;
  auto &block_offset_0 = this->block_offset_0_;
  auto &omega_odd = this->omega_odd_;
  auto &seed = this->seed_;
  auto &omega_even = this->omega_even_;
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(force->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(force->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(force->nrOfGhostLayers()))
  float *RESTRICT const _data_force = force->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs->nrOfGhostLayers()))
  float *RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)pdfs->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()))
  float *RESTRICT _data_pdfs_tmp = pdfs_tmp->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs_tmp->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)pdfs_tmp->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 2))
  const int64_t _size_force_0 = int64_t(int64_c(ci.xSize()) + 2);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 2))
  const int64_t _size_force_1 = int64_t(int64_c(ci.ySize()) + 2);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 2))
  const int64_t _size_force_2 = int64_t(int64_c(ci.zSize()) + 2);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0)
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  const int64_t _stride_pdfs_tmp_1 = int64_t(pdfs_tmp->yStride());
  const int64_t _stride_pdfs_tmp_2 = int64_t(pdfs_tmp->zStride());
  const int64_t _stride_pdfs_tmp_3 = int64_t(1 * int64_t(pdfs_tmp->fStride()));
  internal_482af0085ade5d86c474548ac4b32084::streamcollidesweepthermalizedsingleprecisionavx_streamcollidesweepthermalizedsingleprecisionavx(_data_force, _data_pdfs, _data_pdfs_tmp, _size_force_0, _size_force_1, _size_force_2, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_pdfs_tmp_3, block_offset_0, block_offset_1, block_offset_2, kT, omega_bulk, omega_even, omega_odd, omega_shear, seed, time_step);
  pdfs->swapDataPointers(pdfs_tmp);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif
