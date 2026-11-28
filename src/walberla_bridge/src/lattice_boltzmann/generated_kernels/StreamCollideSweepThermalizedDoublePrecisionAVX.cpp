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
//! \\file StreamCollideSweepThermalizedDoublePrecisionAVX.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 272d4a09ec35da50685afc9586645e1b9984b423

#include <cmath>

#include "StreamCollideSweepThermalizedDoublePrecisionAVX.h"
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

namespace internal_b312c2da0affb0ff73dd81300a784e62 {
static FUNC_PREFIX void streamcollidesweepthermalizeddoubleprecisionavx_streamcollidesweepthermalizeddoubleprecisionavx(double *RESTRICT const _data_force, double *RESTRICT const _data_pdfs, double *RESTRICT _data_pdfs_tmp, int64_t const _size_force_0, int64_t const _size_force_1, int64_t const _size_force_2, int64_t const _stride_force_1, int64_t const _stride_force_2, int64_t const _stride_force_3, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, int64_t const _stride_pdfs_tmp_1, int64_t const _stride_pdfs_tmp_2, int64_t const _stride_pdfs_tmp_3, uint32_t block_offset_0, uint32_t block_offset_1, uint32_t block_offset_2, double kT, double omega_bulk, double omega_even, double omega_odd, double omega_shear, uint32_t seed, uint32_t time_step) {
#pragma omp parallel
  {
    const double xi_20 = omega_bulk * 0.5;
    const double xi_47 = omega_shear * 0.041666666666666664;
    const double xi_51 = omega_bulk * 0.041666666666666664;
    const double xi_62 = omega_shear * 0.125;
    const double xi_97 = 3.7416573867739413;
    const double xi_100 = 5.4772255750516612;
    const double xi_104 = 2.4494897427831779;
    const double xi_107 = 8.3666002653407556;
    const double xi_150 = omega_odd * 0.25;
    const double xi_160 = omega_odd * 0.083333333333333329;
    const double xi_173 = 1.7320508075688772;
    const double xi_217 = omega_shear * 0.25;
    const double xi_223 = omega_odd * 0.041666666666666664;
    const double xi_226 = omega_odd * 0.125;
    const double rr_0 = 0.0;
    const double xi_45 = rr_0 * 0.041666666666666664;
#pragma omp for schedule(static)
    for (int64_t ctr_2 = 1; ctr_2 < _size_force_2 - 1; ctr_2 += 1) {
      for (int64_t ctr_1 = 1; ctr_1 < _size_force_1 - 1; ctr_1 += 1) {
        {
          for (int64_t ctr_0 = 1; ctr_0 < (int64_t)((_size_force_0 - 2) / (4)) * (4) + 1; ctr_0 += 4) {

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
            const __m256d xi_2 = _mm256_add_pd(_mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3 + ctr_0 - 1]), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3 + ctr_0 - 1]));
            const __m256d xi_3 = _mm256_add_pd(_mm256_add_pd(xi_2, _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3 + ctr_0 - 1])), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1]));
            const __m256d xi_4 = _mm256_add_pd(_mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3 + ctr_0]), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3 + ctr_0 + 1]));
            const __m256d xi_5 = _mm256_add_pd(_mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0]), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0]));
            const __m256d xi_6 = _mm256_add_pd(_mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3 + ctr_0]), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0]));
            const __m256d xi_7 = _mm256_add_pd(_mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3 + ctr_0 + 1]), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1]));
            const __m256d xi_8 = _mm256_add_pd(xi_7, _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1]));
            const __m256d xi_9 = _mm256_add_pd(_mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0]), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0]));
            const __m256d xi_11 = _mm256_add_pd(_mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3 + ctr_0 + 1]), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1]));
            const __m256d xi_12 = _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3 + ctr_0 - 1])), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3 + ctr_0 - 1]));
            const __m256d xi_13 = _mm256_add_pd(_mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1]), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0]));
            const __m256d xi_14 = _mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0]));
            const __m256d xi_15 = _mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3 + ctr_0 - 1]));
            const __m256d xi_16 = _mm256_add_pd(_mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1]), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0]));
            const __m256d xi_21 = _mm256_mul_pd(_mm256_set_pd(0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666666), _mm256_loadu_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0]));
            const __m256d xi_22 = _mm256_mul_pd(_mm256_set_pd(0.083333333333333329, 0.083333333333333329, 0.083333333333333329, 0.083333333333333329), _mm256_loadu_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0]));
            const __m256d xi_33 = _mm256_mul_pd(_mm256_set_pd(0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666666), _mm256_load_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0]));
            const __m256d xi_34 = _mm256_mul_pd(_mm256_set_pd(0.083333333333333329, 0.083333333333333329, 0.083333333333333329, 0.083333333333333329), _mm256_load_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0]));
            const __m256d xi_39 = _mm256_mul_pd(_mm256_set_pd(0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666666), _mm256_loadu_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0]));
            const __m256d xi_40 = _mm256_mul_pd(_mm256_set_pd(0.083333333333333329, 0.083333333333333329, 0.083333333333333329, 0.083333333333333329), _mm256_loadu_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0]));
            const __m256d xi_58 = _mm256_mul_pd(_mm256_set_pd(0.25, 0.25, 0.25, 0.25), _mm256_loadu_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0]));
            const __m256d xi_63 = _mm256_mul_pd(_mm256_set_pd(xi_62, xi_62, xi_62, xi_62), _mm256_loadu_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0]));
            const __m256d xi_98 = _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_6_0);
            const __m256d xi_101 = _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_7_0);
            const __m256d xi_103 = _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_2_1);
            const __m256d xi_108 = _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_6_1);
            const __m256d xi_112 = _mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_load_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + ctr_0]));
            const __m256d xi_118 = _mm256_add_pd(xi_11, xi_3);
            const __m256d xi_121 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(2.0, 2.0, 2.0, 2.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0])), _mm256_mul_pd(_mm256_set_pd(2.0, 2.0, 2.0, 2.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_pd(_mm256_set_pd(2.0, 2.0, 2.0, 2.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_pd(_mm256_set_pd(2.0, 2.0, 2.0, 2.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_pd(_mm256_set_pd(5.0, 5.0, 5.0, 5.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3 + ctr_0 + 1]))), _mm256_mul_pd(_mm256_set_pd(5.0, 5.0, 5.0, 5.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3 + ctr_0 - 1])));
            const __m256d xi_124 = _mm256_mul_pd(_mm256_set_pd(2.0, 2.0, 2.0, 2.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1]));
            const __m256d xi_125 = _mm256_mul_pd(_mm256_set_pd(2.0, 2.0, 2.0, 2.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3 + ctr_0 - 1]));
            const __m256d xi_126 = _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(2.0, 2.0, 2.0, 2.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1])), _mm256_mul_pd(_mm256_set_pd(2.0, 2.0, 2.0, 2.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1])));
            const __m256d xi_129 = _mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3 + ctr_0]));
            const __m256d xi_133 = _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_0_1);
            const __m256d xi_137 = _mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1]));
            const __m256d xi_138 = _mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1]));
            const __m256d xi_139 = _mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3 + ctr_0]));
            const __m256d xi_140 = _mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0]));
            const __m256d xi_145 = _mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0]));
            const __m256d xi_146 = _mm256_add_pd(xi_14, _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0]));
            const __m256d xi_147 = _mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0]));
            const __m256d xi_148 = _mm256_add_pd(xi_147, _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3 + ctr_0]));
            const __m256d xi_149 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_145, xi_146), xi_148), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0]));
            const __m256d xi_151 = _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_5_1);
            const __m256d xi_155 = _mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0]));
            const __m256d xi_156 = _mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0]));
            const __m256d xi_157 = _mm256_mul_pd(_mm256_set_pd(2.0, 2.0, 2.0, 2.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3 + ctr_0 - 1]));
            const __m256d xi_158 = _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(2.0, 2.0, 2.0, 2.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3 + ctr_0 - 1])), _mm256_mul_pd(_mm256_set_pd(-2.0, -2.0, -2.0, -2.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3 + ctr_0 + 1])));
            const __m256d xi_159 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(2.0, 2.0, 2.0, 2.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1])), _mm256_mul_pd(xi_157, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_148), xi_155), xi_156), xi_158), xi_5);
            const __m256d xi_161 = _mm256_mul_pd(xi_159, _mm256_set_pd(xi_160, xi_160, xi_160, xi_160));
            const __m256d xi_162 = _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_3_0);
            const __m256d xi_174 = _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_0_0);
            const __m256d xi_185 = _mm256_add_pd(xi_15, _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1]));
            const __m256d xi_186 = _mm256_add_pd(xi_138, xi_185);
            const __m256d xi_187 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3 + ctr_0 + 1])), xi_186), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3 + ctr_0 - 1])), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1]));
            const __m256d xi_188 = _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_4_1);
            const __m256d xi_189 = _mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3 + ctr_0 - 1]));
            const __m256d xi_190 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(2.0, 2.0, 2.0, 2.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1])), _mm256_mul_pd(xi_157, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_158, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_186, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_189, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_7, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)));
            const __m256d xi_191 = _mm256_mul_pd(xi_190, _mm256_set_pd(xi_160, xi_160, xi_160, xi_160));
            const __m256d xi_192 = _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_4_0);
            const __m256d xi_198 = _mm256_add_pd(_mm256_add_pd(xi_145, xi_155), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0]));
            const __m256d xi_199 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_139, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_198, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0])));
            const __m256d xi_200 = _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_5_0);
            const __m256d xi_201 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_124, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_125, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_126), xi_140), xi_198), xi_6);
            const __m256d xi_202 = _mm256_mul_pd(xi_201, _mm256_set_pd(xi_160, xi_160, xi_160, xi_160));
            const __m256d xi_203 = _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_3_1);
            const __m256d xi_224 = _mm256_mul_pd(xi_159, _mm256_set_pd(xi_223, xi_223, xi_223, xi_223));
            const __m256d xi_227 = _mm256_mul_pd(xi_149, _mm256_set_pd(xi_226, xi_226, xi_226, xi_226));
            const __m256d xi_233 = _mm256_mul_pd(xi_201, _mm256_set_pd(xi_223, xi_223, xi_223, xi_223));
            const __m256d xi_234 = _mm256_mul_pd(xi_199, _mm256_set_pd(xi_226, xi_226, xi_226, xi_226));
            const __m256d xi_254 = _mm256_mul_pd(xi_187, _mm256_set_pd(xi_226, xi_226, xi_226, xi_226));
            const __m256d xi_255 = _mm256_mul_pd(xi_190, _mm256_set_pd(xi_223, xi_223, xi_223, xi_223));
            const __m256d xi_23 = _mm256_mul_pd(xi_22, _mm256_set_pd(rr_0, rr_0, rr_0, rr_0));
            const __m256d xi_35 = _mm256_mul_pd(xi_34, _mm256_set_pd(rr_0, rr_0, rr_0, rr_0));
            const __m256d xi_41 = _mm256_mul_pd(xi_40, _mm256_set_pd(rr_0, rr_0, rr_0, rr_0));
            const __m256d xi_46 = _mm256_mul_pd(_mm256_set_pd(xi_45, xi_45, xi_45, xi_45), _mm256_load_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0]));
            const __m256d xi_50 = _mm256_mul_pd(_mm256_set_pd(xi_45, xi_45, xi_45, xi_45), _mm256_loadu_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0]));
            const __m256d xi_72 = _mm256_mul_pd(_mm256_set_pd(xi_45, xi_45, xi_45, xi_45), _mm256_loadu_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0]));
            const __m256d vel0Term = _mm256_add_pd(xi_3, _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3 + ctr_0 - 1]));
            const __m256d vel1Term = _mm256_add_pd(xi_4, xi_5);
            const __m256d vel2Term = _mm256_add_pd(xi_6, _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1]));
            const __m256d delta_rho = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(vel0Term, vel1Term), vel2Term), xi_8), xi_9), _mm256_load_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + ctr_0])), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0]));
            const __m256d rho = _mm256_add_pd(_mm256_set_pd(1.0, 1.0, 1.0, 1.0), delta_rho);
            const __m256d xi_95 = _mm256_mul_pd(rho, _mm256_set_pd(kT, kT, kT, kT));
            const __m256d xi_96 = _mm256_sqrt_pd(_mm256_mul_pd(xi_95, _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_set_pd(omega_even, omega_even, omega_even, omega_even)), _mm256_set_pd(1.0, 1.0, 1.0, 1.0)), _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_set_pd(omega_even, omega_even, omega_even, omega_even)), _mm256_set_pd(1.0, 1.0, 1.0, 1.0)))), _mm256_set_pd(1.0, 1.0, 1.0, 1.0))));
            const __m256d xi_99 = _mm256_mul_pd(_mm256_mul_pd(xi_96, xi_98), _mm256_set_pd(xi_97, xi_97, xi_97, xi_97));
            const __m256d xi_102 = _mm256_mul_pd(_mm256_mul_pd(xi_101, xi_96), _mm256_set_pd(xi_100, xi_100, xi_100, xi_100));
            const __m256d xi_105 = _mm256_sqrt_pd(_mm256_mul_pd(xi_95, _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk, omega_bulk)), _mm256_set_pd(1.0, 1.0, 1.0, 1.0)), _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk, omega_bulk)), _mm256_set_pd(1.0, 1.0, 1.0, 1.0)))), _mm256_set_pd(1.0, 1.0, 1.0, 1.0))));
            const __m256d xi_106 = _mm256_mul_pd(_mm256_mul_pd(xi_103, xi_105), _mm256_set_pd(xi_104, xi_104, xi_104, xi_104));
            const __m256d xi_109 = _mm256_mul_pd(_mm256_mul_pd(xi_108, xi_96), _mm256_set_pd(xi_107, xi_107, xi_107, xi_107));
            const __m256d xi_131 = _mm256_mul_pd(xi_99, _mm256_set_pd(0.11904761904761904, 0.11904761904761904, 0.11904761904761904, 0.11904761904761904));
            const __m256d xi_134 = _mm256_sqrt_pd(_mm256_mul_pd(xi_95, _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_set_pd(omega_shear, omega_shear, omega_shear, omega_shear)), _mm256_set_pd(1.0, 1.0, 1.0, 1.0)), _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_set_pd(omega_shear, omega_shear, omega_shear, omega_shear)), _mm256_set_pd(1.0, 1.0, 1.0, 1.0)))), _mm256_set_pd(1.0, 1.0, 1.0, 1.0))));
            const __m256d xi_135 = _mm256_mul_pd(xi_134, _mm256_set_pd(0.5, 0.5, 0.5, 0.5));
            const __m256d xi_136 = _mm256_mul_pd(xi_133, xi_135);
            const __m256d xi_152 = _mm256_sqrt_pd(_mm256_mul_pd(xi_95, _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_mul_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_set_pd(omega_odd, omega_odd, omega_odd, omega_odd)), _mm256_set_pd(1.0, 1.0, 1.0, 1.0)), _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_set_pd(omega_odd, omega_odd, omega_odd, omega_odd)), _mm256_set_pd(1.0, 1.0, 1.0, 1.0)))), _mm256_set_pd(1.0, 1.0, 1.0, 1.0))));
            const __m256d xi_153 = _mm256_mul_pd(xi_152, _mm256_set_pd(1.4142135623730951, 1.4142135623730951, 1.4142135623730951, 1.4142135623730951));
            const __m256d xi_154 = _mm256_mul_pd(xi_153, _mm256_set_pd(0.5, 0.5, 0.5, 0.5));
            const __m256d xi_163 = _mm256_mul_pd(xi_152, _mm256_set_pd(xi_104, xi_104, xi_104, xi_104));
            const __m256d xi_164 = _mm256_mul_pd(xi_163, _mm256_set_pd(0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666666));
            const __m256d xi_165 = _mm256_mul_pd(xi_162, xi_164);
            const __m256d xi_166 = _mm256_add_pd(xi_161, xi_165);
            const __m256d xi_167 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_149, _mm256_set_pd(xi_150, xi_150, xi_150, xi_150)), _mm256_mul_pd(xi_151, xi_154)), xi_166);
            const __m256d xi_169 = _mm256_mul_pd(xi_102, _mm256_set_pd(0.10000000000000001, 0.10000000000000001, 0.10000000000000001, 0.10000000000000001));
            const __m256d xi_175 = _mm256_mul_pd(_mm256_mul_pd(xi_134, xi_174), _mm256_set_pd(xi_173, xi_173, xi_173, xi_173));
            const __m256d xi_176 = _mm256_mul_pd(xi_175, _mm256_set_pd(0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666666));
            const __m256d xi_184 = _mm256_mul_pd(xi_109, _mm256_set_pd(0.071428571428571425, 0.071428571428571425, 0.071428571428571425, 0.071428571428571425));
            const __m256d xi_193 = _mm256_mul_pd(xi_164, xi_192);
            const __m256d xi_194 = _mm256_add_pd(xi_191, xi_193);
            const __m256d xi_195 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_187, _mm256_set_pd(xi_150, xi_150, xi_150, xi_150)), _mm256_mul_pd(xi_154, xi_188)), xi_194);
            const __m256d xi_197 = _mm256_mul_pd(xi_109, _mm256_set_pd(0.042857142857142858, 0.042857142857142858, 0.042857142857142858, 0.042857142857142858));
            const __m256d xi_204 = _mm256_mul_pd(xi_164, xi_203);
            const __m256d xi_205 = _mm256_add_pd(xi_202, xi_204);
            const __m256d xi_206 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_199, _mm256_set_pd(xi_150, xi_150, xi_150, xi_150)), _mm256_mul_pd(xi_154, xi_200)), xi_205);
            const __m256d xi_207 = _mm256_mul_pd(_mm256_mul_pd(xi_133, xi_134), _mm256_set_pd(0.25, 0.25, 0.25, 0.25));
            const __m256d xi_210 = _mm256_mul_pd(xi_99, _mm256_set_pd(0.083333333333333329, 0.083333333333333329, 0.083333333333333329, 0.083333333333333329));
            const __m256d xi_214 = _mm256_add_pd(_mm256_mul_pd(xi_191, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_193, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)));
            const __m256d xi_215 = _mm256_mul_pd(xi_135, _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_1_0));
            const __m256d xi_222 = _mm256_mul_pd(xi_135, _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_2_0));
            const __m256d xi_228 = _mm256_mul_pd(xi_163, _mm256_set_pd(0.083333333333333329, 0.083333333333333329, 0.083333333333333329, 0.083333333333333329));
            const __m256d xi_229 = _mm256_mul_pd(xi_162, xi_228);
            const __m256d xi_230 = _mm256_mul_pd(xi_153, _mm256_set_pd(0.25, 0.25, 0.25, 0.25));
            const __m256d xi_231 = _mm256_mul_pd(xi_151, xi_230);
            const __m256d xi_235 = _mm256_mul_pd(xi_203, xi_228);
            const __m256d xi_236 = _mm256_mul_pd(xi_200, xi_230);
            const __m256d xi_237 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_233, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_235, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_234), xi_236);
            const __m256d xi_239 = _mm256_mul_pd(xi_109, _mm256_set_pd(0.014285714285714285, 0.014285714285714285, 0.014285714285714285, 0.014285714285714285));
            const __m256d xi_241 = _mm256_mul_pd(xi_99, _mm256_set_pd(0.023809523809523808, 0.023809523809523808, 0.023809523809523808, 0.023809523809523808));
            const __m256d xi_244 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_234, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_236, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_233), xi_235);
            const __m256d xi_246 = _mm256_mul_pd(xi_207, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
            const __m256d xi_249 = _mm256_mul_pd(xi_109, _mm256_set_pd(0.035714285714285712, 0.035714285714285712, 0.035714285714285712, 0.035714285714285712));
            const __m256d xi_251 = _mm256_mul_pd(xi_135, _mm256_add_pd(_mm256_set_pd(-0.5, -0.5, -0.5, -0.5), random_1_1));
            const __m256d xi_256 = _mm256_mul_pd(xi_188, xi_230);
            const __m256d xi_257 = _mm256_mul_pd(xi_192, xi_228);
            const __m256d xi_258 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_254, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_256, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_255), xi_257);
            const __m256d xi_260 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_255, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_257, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_254), xi_256);
            const __m256d xi_0 = _mm256_div_pd(_mm256_set_pd(1.0, 1.0, 1.0, 1.0), rho);
            const __m256d xi_10 = _mm256_mul_pd(xi_0, _mm256_set_pd(0.5, 0.5, 0.5, 0.5));
            const __m256d u_0 = _mm256_add_pd(_mm256_mul_pd(xi_0, _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_11, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_8, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), vel0Term)), _mm256_mul_pd(xi_10, _mm256_load_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0])));
            const __m256d xi_17 = _mm256_mul_pd(u_0, _mm256_load_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0]));
            const __m256d xi_28 = _mm256_mul_pd(xi_17, _mm256_set_pd(0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666666));
            const __m256d xi_29 = _mm256_mul_pd(xi_28, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
            const __m256d xi_30 = _mm256_mul_pd(xi_17, _mm256_set_pd(0.083333333333333329, 0.083333333333333329, 0.083333333333333329, 0.083333333333333329));
            const __m256d xi_31 = _mm256_add_pd(_mm256_mul_pd(xi_30, _mm256_set_pd(omega_shear, omega_shear, omega_shear, omega_shear)), xi_29);
            const __m256d xi_48 = _mm256_add_pd(_mm256_mul_pd(xi_17, _mm256_set_pd(xi_47, xi_47, xi_47, xi_47)), xi_29);
            const __m256d xi_49 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_46, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_34), xi_48);
            const __m256d xi_52 = _mm256_mul_pd(xi_17, _mm256_set_pd(xi_51, xi_51, xi_51, xi_51));
            const __m256d xi_59 = _mm256_mul_pd(u_0, xi_58);
            const __m256d xi_64 = _mm256_mul_pd(u_0, xi_63);
            const __m256d xi_68 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_34, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_46), xi_48);
            const __m256d xi_75 = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(u_0, _mm256_set_pd(-0.083333333333333329, -0.083333333333333329, -0.083333333333333329, -0.083333333333333329)), _mm256_set_pd(omega_shear, omega_shear, omega_shear, omega_shear)), _mm256_load_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0]));
            const __m256d xi_85 = _mm256_mul_pd(u_0, _mm256_loadu_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0]));
            const __m256d xi_86 = _mm256_mul_pd(xi_85, _mm256_set_pd(0.25, 0.25, 0.25, 0.25));
            const __m256d xi_89 = _mm256_mul_pd(xi_85, _mm256_set_pd(xi_62, xi_62, xi_62, xi_62));
            const __m256d xi_111 = _mm256_mul_pd(u_0, u_0);
            const __m256d u_1 = _mm256_add_pd(_mm256_mul_pd(xi_0, _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_12, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_13, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_9, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), vel1Term)), _mm256_mul_pd(xi_10, _mm256_loadu_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0])));
            const __m256d xi_18 = _mm256_mul_pd(u_1, _mm256_loadu_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0]));
            const __m256d xi_26 = _mm256_mul_pd(xi_18, _mm256_set_pd(0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666666));
            const __m256d xi_36 = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(u_1, _mm256_set_pd(-0.083333333333333329, -0.083333333333333329, -0.083333333333333329, -0.083333333333333329)), _mm256_set_pd(omega_shear, omega_shear, omega_shear, omega_shear)), _mm256_loadu_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0]));
            const __m256d xi_42 = _mm256_mul_pd(xi_26, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
            const __m256d xi_43 = _mm256_mul_pd(xi_18, _mm256_set_pd(0.083333333333333329, 0.083333333333333329, 0.083333333333333329, 0.083333333333333329));
            const __m256d xi_53 = _mm256_mul_pd(xi_18, _mm256_set_pd(xi_51, xi_51, xi_51, xi_51));
            const __m256d xi_60 = _mm256_mul_pd(u_1, _mm256_set_pd(0.25, 0.25, 0.25, 0.25));
            const __m256d xi_61 = _mm256_mul_pd(xi_60, _mm256_load_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0]));
            const __m256d xi_65 = _mm256_mul_pd(u_1, _mm256_set_pd(xi_62, xi_62, xi_62, xi_62));
            const __m256d xi_66 = _mm256_mul_pd(xi_65, _mm256_load_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0]));
            const __m256d xi_67 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_64, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_66, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_59), xi_61);
            const __m256d xi_69 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_59, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_61, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_64), xi_66);
            const __m256d xi_77 = _mm256_mul_pd(xi_60, _mm256_loadu_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0]));
            const __m256d xi_79 = _mm256_mul_pd(xi_65, _mm256_loadu_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0]));
            const __m256d xi_110 = _mm256_mul_pd(rho, _mm256_mul_pd(u_1, u_1));
            const __m256d xi_117 = _mm256_mul_pd(xi_110, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
            const __m256d xi_216 = _mm256_mul_pd(rho, u_1);
            const __m256d xi_218 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1])), _mm256_mul_pd(u_0, xi_216)), xi_12), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3 + ctr_0 + 1])), _mm256_set_pd(xi_217, xi_217, xi_217, xi_217));
            const __m256d xi_219 = _mm256_add_pd(_mm256_mul_pd(xi_215, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_218, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)));
            const __m256d xi_220 = _mm256_add_pd(xi_215, xi_218);
            const __m256d u_2 = _mm256_add_pd(_mm256_mul_pd(xi_0, _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_14, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_15, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_16, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1]))), _mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0]))), vel2Term)), _mm256_mul_pd(xi_10, _mm256_loadu_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0])));
            const __m256d xi_19 = _mm256_mul_pd(u_2, _mm256_loadu_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0]));
            const __m256d xi_24 = _mm256_mul_pd(xi_19, _mm256_set_pd(0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666666));
            const __m256d xi_25 = _mm256_mul_pd(xi_24, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
            const __m256d xi_27 = _mm256_mul_pd(xi_19, _mm256_set_pd(0.083333333333333329, 0.083333333333333329, 0.083333333333333329, 0.083333333333333329));
            const __m256d xi_32 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_18, _mm256_set_pd(0.33333333333333331, 0.33333333333333331, 0.33333333333333331, 0.33333333333333331)), _mm256_mul_pd(xi_27, _mm256_set_pd(omega_shear, omega_shear, omega_shear, omega_shear))), _mm256_mul_pd(_mm256_mul_pd(xi_26, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_set_pd(omega_shear, omega_shear, omega_shear, omega_shear))), xi_25), xi_31);
            const __m256d xi_37 = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(u_2, _mm256_set_pd(-0.083333333333333329, -0.083333333333333329, -0.083333333333333329, -0.083333333333333329)), _mm256_set_pd(omega_shear, omega_shear, omega_shear, omega_shear)), _mm256_loadu_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0]));
            const __m256d xi_38 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_28, _mm256_set_pd(omega_shear, omega_shear, omega_shear, omega_shear)), _mm256_mul_pd(_mm256_mul_pd(u_0, _mm256_set_pd(-0.33333333333333331, -0.33333333333333331, -0.33333333333333331, -0.33333333333333331)), _mm256_load_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0]))), xi_24), xi_26), xi_36), xi_37);
            const __m256d xi_44 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_19, _mm256_set_pd(0.33333333333333331, 0.33333333333333331, 0.33333333333333331, 0.33333333333333331)), _mm256_mul_pd(xi_43, _mm256_set_pd(omega_shear, omega_shear, omega_shear, omega_shear))), _mm256_mul_pd(_mm256_mul_pd(xi_24, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_set_pd(omega_shear, omega_shear, omega_shear, omega_shear))), xi_31), xi_42);
            const __m256d xi_54 = _mm256_mul_pd(xi_19, _mm256_set_pd(xi_51, xi_51, xi_51, xi_51));
            const __m256d xi_55 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_18, _mm256_set_pd(xi_47, xi_47, xi_47, xi_47)), xi_42), xi_52), xi_53), xi_54);
            const __m256d xi_56 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_22, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_50), xi_55);
            const __m256d xi_57 = _mm256_add_pd(_mm256_add_pd(xi_27, xi_37), xi_56);
            const __m256d xi_70 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_50, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_22), xi_55);
            const __m256d xi_71 = _mm256_add_pd(_mm256_add_pd(xi_27, xi_37), xi_70);
            const __m256d xi_73 = _mm256_add_pd(_mm256_mul_pd(xi_19, _mm256_set_pd(xi_47, xi_47, xi_47, xi_47)), xi_25);
            const __m256d xi_74 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_40, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_72), xi_73);
            const __m256d xi_76 = _mm256_add_pd(_mm256_add_pd(xi_30, xi_56), xi_75);
            const __m256d xi_78 = _mm256_mul_pd(u_2, xi_58);
            const __m256d xi_80 = _mm256_mul_pd(u_2, xi_63);
            const __m256d xi_81 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_77, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_78, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_79), xi_80);
            const __m256d xi_82 = _mm256_add_pd(_mm256_add_pd(xi_30, xi_70), xi_75);
            const __m256d xi_83 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_79, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_80, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_77), xi_78);
            const __m256d xi_84 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_36, xi_43), xi_52), xi_53), xi_54), xi_74);
            const __m256d xi_87 = _mm256_mul_pd(u_2, _mm256_load_pd(&_data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0]));
            const __m256d xi_88 = _mm256_mul_pd(xi_87, _mm256_set_pd(0.25, 0.25, 0.25, 0.25));
            const __m256d xi_90 = _mm256_mul_pd(xi_87, _mm256_set_pd(xi_62, xi_62, xi_62, xi_62));
            const __m256d xi_91 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_89, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_90, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_86), xi_88);
            const __m256d xi_92 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_86, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_88, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_89), xi_90);
            const __m256d xi_93 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_72, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_40), xi_73);
            const __m256d xi_94 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_36, xi_43), xi_52), xi_53), xi_54), xi_93);
            const __m256d xi_113 = _mm256_mul_pd(rho, _mm256_mul_pd(u_2, u_2));
            const __m256d xi_114 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(3.0, 3.0, 3.0, 3.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3 + ctr_0])), _mm256_mul_pd(_mm256_set_pd(3.0, 3.0, 3.0, 3.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_pd(xi_113, _mm256_set_pd(0.66666666666666663, 0.66666666666666663, 0.66666666666666663, 0.66666666666666663))), xi_112);
            const __m256d xi_115 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(3.0, 3.0, 3.0, 3.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3 + ctr_0])), _mm256_mul_pd(_mm256_set_pd(3.0, 3.0, 3.0, 3.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_pd(xi_110, _mm256_set_pd(0.66666666666666663, 0.66666666666666663, 0.66666666666666663, 0.66666666666666663))), _mm256_mul_pd(_mm256_set_pd(-3.0, -3.0, -3.0, -3.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_pd(_mm256_set_pd(-3.0, -3.0, -3.0, -3.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_pd(_mm256_set_pd(-3.0, -3.0, -3.0, -3.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_pd(_mm256_set_pd(-3.0, -3.0, -3.0, -3.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_pd(_mm256_mul_pd(rho, xi_111), _mm256_set_pd(1.6666666666666667, 1.6666666666666667, 1.6666666666666667, 1.6666666666666667))), xi_114);
            const __m256d xi_116 = _mm256_mul_pd(xi_115, _mm256_set_pd(omega_even, omega_even, omega_even, omega_even));
            const __m256d xi_119 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_112, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_117, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_118, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_13, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_16, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_5, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(rho, xi_111)), xi_113);
            const __m256d xi_120 = _mm256_mul_pd(xi_119, _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk, omega_bulk));
            const __m256d xi_122 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_110, _mm256_set_pd(2.3333333333333335, 2.3333333333333335, 2.3333333333333335, 2.3333333333333335)), _mm256_mul_pd(_mm256_set_pd(-2.0, -2.0, -2.0, -2.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3 + ctr_0]))), _mm256_mul_pd(_mm256_set_pd(-2.0, -2.0, -2.0, -2.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_pd(_mm256_set_pd(-5.0, -5.0, -5.0, -5.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1]))), _mm256_mul_pd(_mm256_set_pd(-5.0, -5.0, -5.0, -5.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1]))), _mm256_mul_pd(_mm256_set_pd(-5.0, -5.0, -5.0, -5.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3 + ctr_0 - 1]))), _mm256_mul_pd(_mm256_set_pd(-5.0, -5.0, -5.0, -5.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1]))), xi_114), xi_121);
            const __m256d xi_123 = _mm256_mul_pd(xi_122, _mm256_set_pd(omega_even, omega_even, omega_even, omega_even));
            const __m256d xi_127 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_113, _mm256_set_pd(3.0, 3.0, 3.0, 3.0)), _mm256_mul_pd(_mm256_set_pd(5.0, 5.0, 5.0, 5.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3 + ctr_0]))), _mm256_mul_pd(_mm256_set_pd(5.0, 5.0, 5.0, 5.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_pd(_mm256_set_pd(-4.0, -4.0, -4.0, -4.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_pd(_mm256_set_pd(-4.0, -4.0, -4.0, -4.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_pd(_mm256_set_pd(-7.0, -7.0, -7.0, -7.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3 + ctr_0 + 1]))), _mm256_mul_pd(_mm256_set_pd(-7.0, -7.0, -7.0, -7.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1]))), _mm256_mul_pd(_mm256_set_pd(-7.0, -7.0, -7.0, -7.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3 + ctr_0 - 1]))), _mm256_mul_pd(_mm256_set_pd(-7.0, -7.0, -7.0, -7.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3 + ctr_0 - 1]))), xi_112), xi_121), xi_124), xi_125), xi_126);
            const __m256d xi_128 = _mm256_mul_pd(xi_127, _mm256_set_pd(omega_even, omega_even, omega_even, omega_even));
            const __m256d xi_130 = _mm256_mul_pd(xi_128, _mm256_set_pd(0.01984126984126984, 0.01984126984126984, 0.01984126984126984, 0.01984126984126984));
            const __m256d xi_132 = _mm256_add_pd(xi_130, xi_131);
            const __m256d xi_141 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_113, xi_139), xi_140), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1]));
            const __m256d xi_142 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_117, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_137, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_138, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_141, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_15, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_2, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_4, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0]))), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1])), _mm256_set_pd(omega_shear, omega_shear, omega_shear, omega_shear));
            const __m256d xi_143 = _mm256_mul_pd(xi_142, _mm256_set_pd(0.125, 0.125, 0.125, 0.125));
            const __m256d xi_144 = _mm256_add_pd(_mm256_mul_pd(xi_136, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_143, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)));
            const __m256d xi_168 = _mm256_mul_pd(xi_116, _mm256_set_pd(0.050000000000000003, 0.050000000000000003, 0.050000000000000003, 0.050000000000000003));
            const __m256d xi_170 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(2.0, 2.0, 2.0, 2.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0])), _mm256_mul_pd(_mm256_set_pd(2.0, 2.0, 2.0, 2.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_pd(_mm256_set_pd(2.0, 2.0, 2.0, 2.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_pd(_mm256_set_pd(2.0, 2.0, 2.0, 2.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0]))), _mm256_mul_pd(xi_110, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_118, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_129, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_141, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_147, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(_mm256_set_pd(-1.0, -1.0, -1.0, -1.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1]))), _mm256_mul_pd(_mm256_set_pd(-2.0, -2.0, -2.0, -2.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3 + ctr_0 + 1]))), _mm256_mul_pd(_mm256_set_pd(-2.0, -2.0, -2.0, -2.0), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3 + ctr_0 - 1]))), _mm256_mul_pd(_mm256_mul_pd(rho, xi_111), _mm256_set_pd(2.0, 2.0, 2.0, 2.0)));
            const __m256d xi_171 = _mm256_mul_pd(xi_170, _mm256_set_pd(omega_shear, omega_shear, omega_shear, omega_shear));
            const __m256d xi_172 = _mm256_mul_pd(xi_171, _mm256_set_pd(0.041666666666666664, 0.041666666666666664, 0.041666666666666664, 0.041666666666666664));
            const __m256d xi_177 = _mm256_add_pd(xi_172, xi_176);
            const __m256d xi_178 = _mm256_add_pd(_mm256_add_pd(xi_168, xi_169), xi_177);
            const __m256d xi_179 = _mm256_add_pd(_mm256_mul_pd(xi_130, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_131, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)));
            const __m256d xi_180 = _mm256_add_pd(xi_136, xi_143);
            const __m256d xi_181 = _mm256_add_pd(_mm256_mul_pd(xi_172, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_176, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)));
            const __m256d xi_182 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_168, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_169, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_181);
            const __m256d xi_183 = _mm256_mul_pd(xi_123, _mm256_set_pd(0.035714285714285712, 0.035714285714285712, 0.035714285714285712, 0.035714285714285712));
            const __m256d xi_196 = _mm256_mul_pd(xi_123, _mm256_set_pd(0.021428571428571429, 0.021428571428571429, 0.021428571428571429, 0.021428571428571429));
            const __m256d xi_208 = _mm256_mul_pd(xi_142, _mm256_set_pd(0.0625, 0.0625, 0.0625, 0.0625));
            const __m256d xi_209 = _mm256_mul_pd(xi_128, _mm256_set_pd(0.013888888888888888, 0.013888888888888888, 0.013888888888888888, 0.013888888888888888));
            const __m256d xi_211 = _mm256_add_pd(_mm256_mul_pd(xi_120, _mm256_set_pd(0.041666666666666664, 0.041666666666666664, 0.041666666666666664, 0.041666666666666664)), _mm256_mul_pd(xi_106, _mm256_set_pd(0.083333333333333329, 0.083333333333333329, 0.083333333333333329, 0.083333333333333329)));
            const __m256d xi_212 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_171, _mm256_set_pd(0.020833333333333332, 0.020833333333333332, 0.020833333333333332, 0.020833333333333332)), _mm256_mul_pd(xi_175, _mm256_set_pd(0.083333333333333329, 0.083333333333333329, 0.083333333333333329, 0.083333333333333329))), xi_211);
            const __m256d xi_213 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_166, xi_207), xi_208), xi_209), xi_210), xi_212);
            const __m256d xi_221 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_161, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_165, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_207), xi_208), xi_209), xi_210), xi_212);
            const __m256d xi_225 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(u_2, xi_216), xi_146), xi_155), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0])), _mm256_set_pd(xi_217, xi_217, xi_217, xi_217));
            const __m256d xi_232 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_224, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_229, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_222), xi_225), xi_227), xi_231);
            const __m256d xi_238 = _mm256_mul_pd(xi_123, _mm256_set_pd(0.0071428571428571426, 0.0071428571428571426, 0.0071428571428571426, 0.0071428571428571426));
            const __m256d xi_240 = _mm256_mul_pd(xi_128, _mm256_set_pd(0.003968253968253968, 0.003968253968253968, 0.003968253968253968, 0.003968253968253968));
            const __m256d xi_242 = _mm256_add_pd(_mm256_mul_pd(xi_240, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_241, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)));
            const __m256d xi_243 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_116, _mm256_set_pd(0.025000000000000001, 0.025000000000000001, 0.025000000000000001, 0.025000000000000001)), _mm256_mul_pd(xi_102, _mm256_set_pd(0.050000000000000003, 0.050000000000000003, 0.050000000000000003, 0.050000000000000003))), _mm256_mul_pd(xi_238, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_239, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_181), xi_211), xi_242);
            const __m256d xi_245 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_115, _mm256_set_pd(-0.025000000000000001, -0.025000000000000001, -0.025000000000000001, -0.025000000000000001)), _mm256_set_pd(omega_even, omega_even, omega_even, omega_even)), _mm256_mul_pd(_mm256_mul_pd(xi_119, _mm256_set_pd(-0.041666666666666664, -0.041666666666666664, -0.041666666666666664, -0.041666666666666664)), _mm256_set_pd(omega_bulk, omega_bulk, omega_bulk, omega_bulk))), _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(xi_101, xi_96), _mm256_set_pd(-0.050000000000000003, -0.050000000000000003, -0.050000000000000003, -0.050000000000000003)), _mm256_set_pd(xi_100, xi_100, xi_100, xi_100))), _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(xi_103, xi_105), _mm256_set_pd(-0.083333333333333329, -0.083333333333333329, -0.083333333333333329, -0.083333333333333329)), _mm256_set_pd(xi_104, xi_104, xi_104, xi_104))), xi_177), xi_238), xi_239), xi_240), xi_241);
            const __m256d xi_247 = _mm256_mul_pd(xi_208, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0));
            const __m256d xi_248 = _mm256_mul_pd(xi_123, _mm256_set_pd(0.017857142857142856, 0.017857142857142856, 0.017857142857142856, 0.017857142857142856));
            const __m256d xi_250 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_205, xi_212), xi_242), xi_246), xi_247), xi_248), xi_249);
            const __m256d xi_252 = _mm256_mul_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(rho, u_0), u_2), xi_137), xi_185), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1])), _mm256_set_pd(xi_217, xi_217, xi_217, xi_217));
            const __m256d xi_253 = _mm256_add_pd(_mm256_mul_pd(xi_251, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_252, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)));
            const __m256d xi_259 = _mm256_add_pd(xi_251, xi_252);
            const __m256d xi_261 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_227, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_231, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_222), xi_224), xi_225), xi_229);
            const __m256d xi_262 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_202, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_204, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_212), xi_242), xi_246), xi_247), xi_248), xi_249);
            const __m256d forceTerm_0 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_17, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_18, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_19, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_17, _mm256_set_pd(xi_20, xi_20, xi_20, xi_20))), _mm256_mul_pd(xi_18, _mm256_set_pd(xi_20, xi_20, xi_20, xi_20))), _mm256_mul_pd(xi_19, _mm256_set_pd(xi_20, xi_20, xi_20, xi_20)));
            const __m256d forceTerm_1 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_23, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_21), xi_32);
            const __m256d forceTerm_2 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_21, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_23), xi_32);
            const __m256d forceTerm_3 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_33, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_38, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_35);
            const __m256d forceTerm_4 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_35, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_38, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), xi_33);
            const __m256d forceTerm_5 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_41, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_39), xi_44);
            const __m256d forceTerm_6 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_39, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), xi_41), xi_44);
            const __m256d forceTerm_7 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_49, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_57, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_67, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)));
            const __m256d forceTerm_8 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_57, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_68, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_69, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)));
            const __m256d forceTerm_9 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_49, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_69, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_71, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)));
            const __m256d forceTerm_10 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_67, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_68, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_71, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)));
            const __m256d forceTerm_11 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_74, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_76, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_81, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)));
            const __m256d forceTerm_12 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_74, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_82, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_83, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)));
            const __m256d forceTerm_13 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_49, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_84, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_91, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)));
            const __m256d forceTerm_14 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_68, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_84, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_92, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)));
            const __m256d forceTerm_15 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_76, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_83, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_93, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)));
            const __m256d forceTerm_16 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_81, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_82, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_93, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)));
            const __m256d forceTerm_17 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_49, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_92, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_94, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)));
            const __m256d forceTerm_18 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_68, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_91, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_94, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)));
            _mm256_store_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_128, _mm256_set_pd(0.023809523809523808, 0.023809523809523808, 0.023809523809523808, 0.023809523809523808)), _mm256_mul_pd(xi_99, _mm256_set_pd(0.14285714285714285, 0.14285714285714285, 0.14285714285714285, 0.14285714285714285))), _mm256_mul_pd(xi_123, _mm256_set_pd(0.042857142857142858, 0.042857142857142858, 0.042857142857142858, 0.042857142857142858))), _mm256_mul_pd(xi_109, _mm256_set_pd(0.085714285714285715, 0.085714285714285715, 0.085714285714285715, 0.085714285714285715))), _mm256_mul_pd(xi_116, _mm256_set_pd(0.10000000000000001, 0.10000000000000001, 0.10000000000000001, 0.10000000000000001))), _mm256_mul_pd(xi_102, _mm256_set_pd(0.20000000000000001, 0.20000000000000001, 0.20000000000000001, 0.20000000000000001))), _mm256_mul_pd(xi_120, _mm256_set_pd(-0.5, -0.5, -0.5, -0.5))), _mm256_mul_pd(xi_106, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), forceTerm_0), _mm256_load_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + ctr_0])));
            _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + _stride_pdfs_tmp_3 + ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_129, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_132, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_144, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_167, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_178, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(_mm256_mul_pd(xi_122, _mm256_set_pd(0.014285714285714285, 0.014285714285714285, 0.014285714285714285, 0.014285714285714285)), _mm256_set_pd(omega_even, omega_even, omega_even, omega_even))), _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(xi_108, xi_96), _mm256_set_pd(0.028571428571428571, 0.028571428571428571, 0.028571428571428571, 0.028571428571428571)), _mm256_set_pd(xi_107, xi_107, xi_107, xi_107))), forceTerm_1));
            _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 2 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_123, _mm256_set_pd(0.014285714285714285, 0.014285714285714285, 0.014285714285714285, 0.014285714285714285)), _mm256_mul_pd(xi_109, _mm256_set_pd(0.028571428571428571, 0.028571428571428571, 0.028571428571428571, 0.028571428571428571))), forceTerm_2), xi_167), xi_179), xi_180), xi_182), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0])));
            _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 3 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_171, _mm256_set_pd(0.083333333333333329, 0.083333333333333329, 0.083333333333333329, 0.083333333333333329)), _mm256_mul_pd(xi_175, _mm256_set_pd(0.33333333333333331, 0.33333333333333331, 0.33333333333333331, 0.33333333333333331))), _mm256_mul_pd(xi_183, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_184, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), forceTerm_3), xi_179), xi_195), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3 + ctr_0 + 1])));
            _mm256_store_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 4 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_132, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_183, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_184, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_189, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_195, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(_mm256_mul_pd(xi_170, _mm256_set_pd(0.083333333333333329, 0.083333333333333329, 0.083333333333333329, 0.083333333333333329)), _mm256_set_pd(omega_shear, omega_shear, omega_shear, omega_shear))), _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(xi_134, xi_174), _mm256_set_pd(0.33333333333333331, 0.33333333333333331, 0.33333333333333331, 0.33333333333333331)), _mm256_set_pd(xi_173, xi_173, xi_173, xi_173))), forceTerm_4));
            _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 5 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_139, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_178, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_180, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_196, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_197, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_206, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(_mm256_mul_pd(xi_127, _mm256_set_pd(0.015873015873015872, 0.015873015873015872, 0.015873015873015872, 0.015873015873015872)), _mm256_set_pd(omega_even, omega_even, omega_even, omega_even))), _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(xi_96, xi_98), _mm256_set_pd(0.095238095238095233, 0.095238095238095233, 0.095238095238095233, 0.095238095238095233)), _mm256_set_pd(xi_97, xi_97, xi_97, xi_97))), forceTerm_5));
            _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 6 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_128, _mm256_set_pd(0.015873015873015872, 0.015873015873015872, 0.015873015873015872, 0.015873015873015872)), _mm256_mul_pd(xi_99, _mm256_set_pd(0.095238095238095233, 0.095238095238095233, 0.095238095238095233, 0.095238095238095233))), _mm256_mul_pd(xi_196, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_197, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), forceTerm_6), xi_144), xi_182), xi_206), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0])));
            _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 7 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_7, xi_213), xi_214), xi_219), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3 + ctr_0 + 1])));
            _mm256_store_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 8 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_8, xi_194), xi_213), xi_220), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3 + ctr_0 - 1])));
            _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 9 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_9, xi_214), xi_220), xi_221), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1])));
            _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 10 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_10, xi_194), xi_219), xi_221), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3 + ctr_0 - 1])));
            _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 11 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_11, xi_232), xi_237), xi_243), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0])));
            _mm256_store_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 12 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_156, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_232, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_244, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_245, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), forceTerm_12));
            _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 13 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_13, xi_250), xi_253), xi_258), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1])));
            _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 14 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_14, xi_250), xi_259), xi_260), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3 + ctr_0 - 1])));
            _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 15 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_145, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0)), _mm256_mul_pd(xi_237, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_245, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), _mm256_mul_pd(xi_261, _mm256_set_pd(-1.0, -1.0, -1.0, -1.0))), forceTerm_15));
            _mm256_store_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 16 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_16, xi_243), xi_244), xi_261), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0])));
            _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 17 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_17, xi_258), xi_259), xi_262), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1])));
            _mm256_storeu_pd(&_data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 18 * _stride_pdfs_tmp_3 + ctr_0], _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_18, xi_253), xi_260), xi_262), _mm256_loadu_pd(&_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1])));
          }
          for (int64_t ctr_0 = (int64_t)((_size_force_0 - 2) / (4)) * (4) + 1; ctr_0 < _size_force_0 - 1; ctr_0 += 1) {

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
            const double xi_2 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3 + ctr_0 - 1] + _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3 + ctr_0 - 1];
            const double xi_3 = xi_2 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1] + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3 + ctr_0 - 1];
            const double xi_4 = _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3 + ctr_0 + 1] + _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3 + ctr_0];
            const double xi_5 = _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0] + _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0];
            const double xi_6 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0] + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3 + ctr_0];
            const double xi_7 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3 + ctr_0 + 1] + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1];
            const double xi_8 = xi_7 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1];
            const double xi_9 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0] + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0];
            const double xi_11 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1] + _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3 + ctr_0 + 1];
            const double xi_12 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3 + ctr_0 - 1] - _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3 + ctr_0 - 1];
            const double xi_13 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1] + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0];
            const double xi_14 = -_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0];
            const double xi_15 = -_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3 + ctr_0 - 1];
            const double xi_16 = _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0] + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1];
            const double xi_21 = 0.16666666666666666 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0];
            const double xi_22 = 0.083333333333333329 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0];
            const double xi_33 = 0.16666666666666666 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0];
            const double xi_34 = 0.083333333333333329 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0];
            const double xi_39 = 0.16666666666666666 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0];
            const double xi_40 = 0.083333333333333329 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0];
            const double xi_58 = 0.25 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0];
            const double xi_63 = xi_62 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0];
            const double xi_98 = random_6_0 - 0.5;
            const double xi_101 = random_7_0 - 0.5;
            const double xi_103 = random_2_1 - 0.5;
            const double xi_108 = random_6_1 - 0.5;
            const double xi_112 = -_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + ctr_0];
            const double xi_118 = xi_11 + xi_3;
            const double xi_121 = 2.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0] + 2.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0] + 2.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0] + 2.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0] + 5.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3 + ctr_0 + 1] + 5.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3 + ctr_0 - 1];
            const double xi_124 = 2.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1];
            const double xi_125 = 2.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3 + ctr_0 - 1];
            const double xi_126 = 2.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1] + 2.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1];
            const double xi_129 = -_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3 + ctr_0];
            const double xi_133 = random_0_1 - 0.5;
            const double xi_137 = -_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1];
            const double xi_138 = -_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1];
            const double xi_139 = -_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3 + ctr_0];
            const double xi_140 = -_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0];
            const double xi_145 = -_data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0];
            const double xi_146 = xi_14 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0];
            const double xi_147 = -_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0];
            const double xi_148 = xi_147 + _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3 + ctr_0];
            const double xi_149 = xi_145 + xi_146 + xi_148 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0];
            const double xi_151 = random_5_1 - 0.5;
            const double xi_155 = -_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0];
            const double xi_156 = -_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0];
            const double xi_157 = 2.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3 + ctr_0 - 1];
            const double xi_158 = -2.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3 + ctr_0 + 1] + 2.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3 + ctr_0 - 1];
            const double xi_159 = xi_148 + xi_155 + xi_156 - xi_157 + xi_158 + xi_5 + 2.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1];
            const double xi_161 = xi_159 * xi_160;
            const double xi_162 = random_3_0 - 0.5;
            const double xi_174 = random_0_0 - 0.5;
            const double xi_185 = xi_15 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1];
            const double xi_186 = xi_138 + xi_185;
            const double xi_187 = xi_186 - _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3 + ctr_0 + 1] + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3 + ctr_0 - 1] + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1];
            const double xi_188 = random_4_1 - 0.5;
            const double xi_189 = -_data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3 + ctr_0 - 1];
            const double xi_190 = -xi_157 - xi_158 - xi_186 - xi_189 - xi_7 + 2.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1];
            const double xi_191 = xi_160 * xi_190;
            const double xi_192 = random_4_0 - 0.5;
            const double xi_198 = xi_145 + xi_155 + _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0];
            const double xi_199 = -xi_139 - xi_198 - _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0] - _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0];
            const double xi_200 = random_5_0 - 0.5;
            const double xi_201 = -xi_124 - xi_125 + xi_126 + xi_140 + xi_198 + xi_6;
            const double xi_202 = xi_160 * xi_201;
            const double xi_203 = random_3_1 - 0.5;
            const double xi_224 = xi_159 * xi_223;
            const double xi_227 = xi_149 * xi_226;
            const double xi_233 = xi_201 * xi_223;
            const double xi_234 = xi_199 * xi_226;
            const double xi_254 = xi_187 * xi_226;
            const double xi_255 = xi_190 * xi_223;
            const double xi_23 = rr_0 * xi_22;
            const double xi_35 = rr_0 * xi_34;
            const double xi_41 = rr_0 * xi_40;
            const double xi_46 = xi_45 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0];
            const double xi_50 = xi_45 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0];
            const double xi_72 = xi_45 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0];
            const double vel0Term = xi_3 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3 + ctr_0 - 1];
            const double vel1Term = xi_4 + xi_5;
            const double vel2Term = xi_6 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1];
            const double delta_rho = vel0Term + vel1Term + vel2Term + xi_8 + xi_9 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0] + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + ctr_0];
            const double rho = delta_rho + 1.0;
            const double xi_95 = kT * rho;
            const double xi_96 = pow(xi_95 * (1.0 - ((-omega_even + 1.0) * (-omega_even + 1.0))), 0.5);
            const double xi_99 = xi_96 * xi_97 * xi_98;
            const double xi_102 = xi_100 * xi_101 * xi_96;
            const double xi_105 = pow(xi_95 * (1.0 - ((-omega_bulk + 1.0) * (-omega_bulk + 1.0))), 0.5);
            const double xi_106 = xi_103 * xi_104 * xi_105;
            const double xi_109 = xi_107 * xi_108 * xi_96;
            const double xi_131 = xi_99 * 0.11904761904761904;
            const double xi_134 = pow(xi_95 * (1.0 - ((-omega_shear + 1.0) * (-omega_shear + 1.0))), 0.5);
            const double xi_135 = xi_134 * 0.5;
            const double xi_136 = xi_133 * xi_135;
            const double xi_152 = pow(xi_95 * (1.0 - ((-omega_odd + 1.0) * (-omega_odd + 1.0))), 0.5);
            const double xi_153 = xi_152 * 1.4142135623730951;
            const double xi_154 = xi_153 * 0.5;
            const double xi_163 = xi_104 * xi_152;
            const double xi_164 = xi_163 * 0.16666666666666666;
            const double xi_165 = xi_162 * xi_164;
            const double xi_166 = xi_161 + xi_165;
            const double xi_167 = xi_149 * xi_150 + xi_151 * xi_154 + xi_166;
            const double xi_169 = xi_102 * 0.10000000000000001;
            const double xi_175 = xi_134 * xi_173 * xi_174;
            const double xi_176 = xi_175 * 0.16666666666666666;
            const double xi_184 = xi_109 * 0.071428571428571425;
            const double xi_193 = xi_164 * xi_192;
            const double xi_194 = xi_191 + xi_193;
            const double xi_195 = xi_150 * xi_187 + xi_154 * xi_188 + xi_194;
            const double xi_197 = xi_109 * 0.042857142857142858;
            const double xi_204 = xi_164 * xi_203;
            const double xi_205 = xi_202 + xi_204;
            const double xi_206 = xi_150 * xi_199 + xi_154 * xi_200 + xi_205;
            const double xi_207 = xi_133 * xi_134 * 0.25;
            const double xi_210 = xi_99 * 0.083333333333333329;
            const double xi_214 = -xi_191 - xi_193;
            const double xi_215 = xi_135 * (random_1_0 - 0.5);
            const double xi_222 = xi_135 * (random_2_0 - 0.5);
            const double xi_228 = xi_163 * 0.083333333333333329;
            const double xi_229 = xi_162 * xi_228;
            const double xi_230 = xi_153 * 0.25;
            const double xi_231 = xi_151 * xi_230;
            const double xi_235 = xi_203 * xi_228;
            const double xi_236 = xi_200 * xi_230;
            const double xi_237 = -xi_233 + xi_234 - xi_235 + xi_236;
            const double xi_239 = xi_109 * 0.014285714285714285;
            const double xi_241 = xi_99 * 0.023809523809523808;
            const double xi_244 = xi_233 - xi_234 + xi_235 - xi_236;
            const double xi_246 = -xi_207;
            const double xi_249 = xi_109 * 0.035714285714285712;
            const double xi_251 = xi_135 * (random_1_1 - 0.5);
            const double xi_256 = xi_188 * xi_230;
            const double xi_257 = xi_192 * xi_228;
            const double xi_258 = -xi_254 + xi_255 - xi_256 + xi_257;
            const double xi_260 = xi_254 - xi_255 + xi_256 - xi_257;
            const double xi_0 = ((1.0) / (rho));
            const double xi_10 = xi_0 * 0.5;
            const double u_0 = xi_0 * (vel0Term - xi_11 - xi_8) + xi_10 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0];
            const double xi_17 = u_0 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0];
            const double xi_28 = xi_17 * 0.16666666666666666;
            const double xi_29 = -xi_28;
            const double xi_30 = xi_17 * 0.083333333333333329;
            const double xi_31 = omega_shear * xi_30 + xi_29;
            const double xi_48 = xi_17 * xi_47 + xi_29;
            const double xi_49 = xi_34 - xi_46 + xi_48;
            const double xi_52 = xi_17 * xi_51;
            const double xi_59 = u_0 * xi_58;
            const double xi_64 = u_0 * xi_63;
            const double xi_68 = -xi_34 + xi_46 + xi_48;
            const double xi_75 = omega_shear * u_0 * -0.083333333333333329 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0];
            const double xi_85 = u_0 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0];
            const double xi_86 = xi_85 * 0.25;
            const double xi_89 = xi_62 * xi_85;
            const double xi_111 = (u_0 * u_0);
            const double u_1 = xi_0 * (vel1Term - xi_12 - xi_13 - xi_9) + xi_10 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0];
            const double xi_18 = u_1 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0];
            const double xi_26 = xi_18 * 0.16666666666666666;
            const double xi_36 = omega_shear * u_1 * -0.083333333333333329 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + _stride_force_3 + ctr_0];
            const double xi_42 = -xi_26;
            const double xi_43 = xi_18 * 0.083333333333333329;
            const double xi_53 = xi_18 * xi_51;
            const double xi_60 = u_1 * 0.25;
            const double xi_61 = xi_60 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0];
            const double xi_65 = u_1 * xi_62;
            const double xi_66 = xi_65 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0];
            const double xi_67 = xi_59 + xi_61 - xi_64 - xi_66;
            const double xi_69 = -xi_59 - xi_61 + xi_64 + xi_66;
            const double xi_77 = xi_60 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0];
            const double xi_79 = xi_65 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0];
            const double xi_110 = rho * (u_1 * u_1);
            const double xi_117 = -xi_110;
            const double xi_216 = rho * u_1;
            const double xi_218 = xi_217 * (u_0 * xi_216 + xi_12 - _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1] + _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3 + ctr_0 + 1]);
            const double xi_219 = -xi_215 - xi_218;
            const double xi_220 = xi_215 + xi_218;
            const double u_2 = xi_0 * (vel2Term - xi_14 - xi_15 - xi_16 - _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1] - _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0] - _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0]) + xi_10 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0];
            const double xi_19 = u_2 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0];
            const double xi_24 = xi_19 * 0.16666666666666666;
            const double xi_25 = -xi_24;
            const double xi_27 = xi_19 * 0.083333333333333329;
            const double xi_32 = -omega_shear * xi_26 + omega_shear * xi_27 + xi_18 * 0.33333333333333331 + xi_25 + xi_31;
            const double xi_37 = omega_shear * u_2 * -0.083333333333333329 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + 2 * _stride_force_3 + ctr_0];
            const double xi_38 = omega_shear * xi_28 + u_0 * -0.33333333333333331 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0] + xi_24 + xi_26 + xi_36 + xi_37;
            const double xi_44 = -omega_shear * xi_24 + omega_shear * xi_43 + xi_19 * 0.33333333333333331 + xi_31 + xi_42;
            const double xi_54 = xi_19 * xi_51;
            const double xi_55 = xi_18 * xi_47 + xi_42 + xi_52 + xi_53 + xi_54;
            const double xi_56 = -xi_22 + xi_50 + xi_55;
            const double xi_57 = xi_27 + xi_37 + xi_56;
            const double xi_70 = xi_22 - xi_50 + xi_55;
            const double xi_71 = xi_27 + xi_37 + xi_70;
            const double xi_73 = xi_19 * xi_47 + xi_25;
            const double xi_74 = -xi_40 + xi_72 + xi_73;
            const double xi_76 = xi_30 + xi_56 + xi_75;
            const double xi_78 = u_2 * xi_58;
            const double xi_80 = u_2 * xi_63;
            const double xi_81 = -xi_77 - xi_78 + xi_79 + xi_80;
            const double xi_82 = xi_30 + xi_70 + xi_75;
            const double xi_83 = xi_77 + xi_78 - xi_79 - xi_80;
            const double xi_84 = xi_36 + xi_43 + xi_52 + xi_53 + xi_54 + xi_74;
            const double xi_87 = u_2 * _data_force[_stride_force_1 * ctr_1 + _stride_force_2 * ctr_2 + ctr_0];
            const double xi_88 = xi_87 * 0.25;
            const double xi_90 = xi_62 * xi_87;
            const double xi_91 = xi_86 + xi_88 - xi_89 - xi_90;
            const double xi_92 = -xi_86 - xi_88 + xi_89 + xi_90;
            const double xi_93 = xi_40 - xi_72 + xi_73;
            const double xi_94 = xi_36 + xi_43 + xi_52 + xi_53 + xi_54 + xi_93;
            const double xi_113 = rho * (u_2 * u_2);
            const double xi_114 = xi_112 + xi_113 * 0.66666666666666663 + 3.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0] + 3.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3 + ctr_0];
            const double xi_115 = rho * xi_111 * 1.6666666666666667 + xi_110 * 0.66666666666666663 + xi_114 - 3.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0] - 3.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0] - 3.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0] - 3.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0] + 3.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0] + 3.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3 + ctr_0];
            const double xi_116 = omega_even * xi_115;
            const double xi_119 = rho * xi_111 - xi_112 + xi_113 - xi_117 - xi_118 - xi_13 - xi_16 - xi_5;
            const double xi_120 = omega_bulk * xi_119;
            const double xi_122 = xi_110 * 2.3333333333333335 + xi_114 + xi_121 - 2.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0] - 2.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3 + ctr_0] - 5.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1] - 5.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1] - 5.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1] - 5.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 14 * _stride_pdfs_3 + ctr_0 - 1];
            const double xi_123 = omega_even * xi_122;
            const double xi_127 = xi_112 + xi_113 * 3.0 + xi_121 + xi_124 + xi_125 + xi_126 - 4.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0] - 4.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 5 * _stride_pdfs_3 + ctr_0] - 7.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3 + ctr_0 - 1] - 7.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1] - 7.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3 + ctr_0 + 1] - 7.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3 + ctr_0 - 1] + 5.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0] + 5.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3 + ctr_0];
            const double xi_128 = omega_even * xi_127;
            const double xi_130 = xi_128 * 0.01984126984126984;
            const double xi_132 = xi_130 + xi_131;
            const double xi_141 = xi_113 + xi_139 + xi_140 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3 + ctr_0 + 1];
            const double xi_142 = omega_shear * (-xi_117 - xi_137 - xi_138 - xi_141 - xi_15 - xi_2 - xi_4 - _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0] + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 13 * _stride_pdfs_3 + ctr_0 + 1]);
            const double xi_143 = xi_142 * 0.125;
            const double xi_144 = -xi_136 - xi_143;
            const double xi_168 = xi_116 * 0.050000000000000003;
            const double xi_170 = rho * xi_111 * 2.0 - xi_110 - xi_118 - xi_129 - xi_141 - xi_147 - 2.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3 + ctr_0 + 1] - 2.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3 + ctr_0 - 1] + 2.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 16 * _stride_pdfs_3 + ctr_0] + 2.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 12 * _stride_pdfs_3 + ctr_0] + 2.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0] + 2.0 * _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 - _stride_pdfs_2 + 11 * _stride_pdfs_3 + ctr_0] - _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 17 * _stride_pdfs_3 + ctr_0 + 1];
            const double xi_171 = omega_shear * xi_170;
            const double xi_172 = xi_171 * 0.041666666666666664;
            const double xi_177 = xi_172 + xi_176;
            const double xi_178 = xi_168 + xi_169 + xi_177;
            const double xi_179 = -xi_130 - xi_131;
            const double xi_180 = xi_136 + xi_143;
            const double xi_181 = -xi_172 - xi_176;
            const double xi_182 = -xi_168 - xi_169 + xi_181;
            const double xi_183 = xi_123 * 0.035714285714285712;
            const double xi_196 = xi_123 * 0.021428571428571429;
            const double xi_208 = xi_142 * 0.0625;
            const double xi_209 = xi_128 * 0.013888888888888888;
            const double xi_211 = xi_106 * 0.083333333333333329 + xi_120 * 0.041666666666666664;
            const double xi_212 = xi_171 * 0.020833333333333332 + xi_175 * 0.083333333333333329 + xi_211;
            const double xi_213 = xi_166 + xi_207 + xi_208 + xi_209 + xi_210 + xi_212;
            const double xi_221 = -xi_161 - xi_165 + xi_207 + xi_208 + xi_209 + xi_210 + xi_212;
            const double xi_225 = xi_217 * (u_2 * xi_216 + xi_146 + xi_155 + _data_pdfs[_stride_pdfs_1 * ctr_1 - _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 15 * _stride_pdfs_3 + ctr_0]);
            const double xi_232 = xi_222 - xi_224 + xi_225 + xi_227 - xi_229 + xi_231;
            const double xi_238 = xi_123 * 0.0071428571428571426;
            const double xi_240 = xi_128 * 0.003968253968253968;
            const double xi_242 = -xi_240 - xi_241;
            const double xi_243 = xi_102 * 0.050000000000000003 + xi_116 * 0.025000000000000001 + xi_181 + xi_211 - xi_238 - xi_239 + xi_242;
            const double xi_245 = omega_bulk * xi_119 * -0.041666666666666664 + omega_even * xi_115 * -0.025000000000000001 + xi_100 * xi_101 * xi_96 * -0.050000000000000003 + xi_103 * xi_104 * xi_105 * -0.083333333333333329 + xi_177 + xi_238 + xi_239 + xi_240 + xi_241;
            const double xi_247 = -xi_208;
            const double xi_248 = xi_123 * 0.017857142857142856;
            const double xi_250 = xi_205 + xi_212 + xi_242 + xi_246 + xi_247 + xi_248 + xi_249;
            const double xi_252 = xi_217 * (rho * u_0 * u_2 + xi_137 + xi_185 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 18 * _stride_pdfs_3 + ctr_0 - 1]);
            const double xi_253 = -xi_251 - xi_252;
            const double xi_259 = xi_251 + xi_252;
            const double xi_261 = xi_222 + xi_224 + xi_225 - xi_227 + xi_229 - xi_231;
            const double xi_262 = -xi_202 - xi_204 + xi_212 + xi_242 + xi_246 + xi_247 + xi_248 + xi_249;
            const double forceTerm_0 = xi_17 * xi_20 - xi_17 + xi_18 * xi_20 - xi_18 + xi_19 * xi_20 - xi_19;
            const double forceTerm_1 = xi_21 - xi_23 + xi_32;
            const double forceTerm_2 = -xi_21 + xi_23 + xi_32;
            const double forceTerm_3 = -xi_33 + xi_35 - xi_38;
            const double forceTerm_4 = xi_33 - xi_35 - xi_38;
            const double forceTerm_5 = xi_39 - xi_41 + xi_44;
            const double forceTerm_6 = -xi_39 + xi_41 + xi_44;
            const double forceTerm_7 = -xi_49 - xi_57 - xi_67;
            const double forceTerm_8 = -xi_57 - xi_68 - xi_69;
            const double forceTerm_9 = -xi_49 - xi_69 - xi_71;
            const double forceTerm_10 = -xi_67 - xi_68 - xi_71;
            const double forceTerm_11 = -xi_74 - xi_76 - xi_81;
            const double forceTerm_12 = -xi_74 - xi_82 - xi_83;
            const double forceTerm_13 = -xi_49 - xi_84 - xi_91;
            const double forceTerm_14 = -xi_68 - xi_84 - xi_92;
            const double forceTerm_15 = -xi_76 - xi_83 - xi_93;
            const double forceTerm_16 = -xi_81 - xi_82 - xi_93;
            const double forceTerm_17 = -xi_49 - xi_92 - xi_94;
            const double forceTerm_18 = -xi_68 - xi_91 - xi_94;
            _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + ctr_0] = forceTerm_0 + xi_102 * 0.20000000000000001 - xi_106 + xi_109 * 0.085714285714285715 + xi_116 * 0.10000000000000001 + xi_120 * -0.5 + xi_123 * 0.042857142857142858 + xi_128 * 0.023809523809523808 + xi_99 * 0.14285714285714285 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + ctr_0];
            _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + _stride_pdfs_tmp_3 + ctr_0] = forceTerm_1 + omega_even * xi_122 * 0.014285714285714285 + xi_107 * xi_108 * xi_96 * 0.028571428571428571 - xi_129 - xi_132 - xi_144 - xi_167 - xi_178;
            _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 2 * _stride_pdfs_tmp_3 + ctr_0] = forceTerm_2 + xi_109 * 0.028571428571428571 + xi_123 * 0.014285714285714285 + xi_167 + xi_179 + xi_180 + xi_182 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_1 + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3 + ctr_0];
            _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 3 * _stride_pdfs_tmp_3 + ctr_0] = forceTerm_3 + xi_171 * 0.083333333333333329 + xi_175 * 0.33333333333333331 + xi_179 - xi_183 - xi_184 + xi_195 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3 + ctr_0 + 1];
            _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 4 * _stride_pdfs_tmp_3 + ctr_0] = forceTerm_4 + omega_shear * xi_170 * 0.083333333333333329 - xi_132 + xi_134 * xi_173 * xi_174 * 0.33333333333333331 - xi_183 - xi_184 - xi_189 - xi_195;
            _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 5 * _stride_pdfs_tmp_3 + ctr_0] = forceTerm_5 + omega_even * xi_127 * 0.015873015873015872 - xi_139 - xi_178 - xi_180 - xi_196 - xi_197 - xi_206 + xi_96 * xi_97 * xi_98 * 0.095238095238095233;
            _data_pdfs_tmp[_stride_pdfs_tmp_1 * ctr_1 + _stride_pdfs_tmp_2 * ctr_2 + 6 * _stride_pdfs_tmp_3 + ctr_0] = forceTerm_6 + xi_128 * 0.015873015873015872 + xi_144 + xi_182 - xi_196 - xi_197 + xi_206 + xi_99 * 0.095238095238095233 + _data_pdfs[_stride_pdfs_1 * ctr_1 + _stride_pdfs_2 * ctr_2 + _stride_pdfs_2 + 6 * _stride_pdfs_3 + ctr_0];
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
} // namespace internal_b312c2da0affb0ff73dd81300a784e62

void StreamCollideSweepThermalizedDoublePrecisionAVX::run(IBlock *block) {
  if (!this->configured_)
    WALBERLA_ABORT("This Sweep contains a configure function that needs to be called manually")

  auto force = block->getData<field::GhostLayerField<double, 3>>(forceID);
  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);
  field::GhostLayerField<double, 19> *pdfs_tmp;
  {
    if (cache_pdfs_.find(block) == cache_pdfs_.end()) {
      pdfs_tmp = pdfs->cloneUninitialized();
      cache_pdfs_[block] = pdfs_tmp;
    } else {
      pdfs_tmp = cache_pdfs_[block];
    }
  }

  auto &seed = this->seed_;
  auto &kT = this->kT_;
  auto &omega_bulk = this->omega_bulk_;
  auto &omega_even = this->omega_even_;
  auto &omega_shear = this->omega_shear_;
  auto &omega_odd = this->omega_odd_;
  auto &block_offset_2 = this->block_offset_2_;
  auto &time_step = this->time_step_;
  auto &block_offset_0 = this->block_offset_0_;
  auto &block_offset_1 = this->block_offset_1_;
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force->nrOfGhostLayers()))
  double *RESTRICT const _data_force = force->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs->nrOfGhostLayers()))
  double *RESTRICT const _data_pdfs = pdfs->dataAt(-1, -1, -1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)pdfs->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_tmp->nrOfGhostLayers()))
  double *RESTRICT _data_pdfs_tmp = pdfs_tmp->dataAt(-1, -1, -1, 0);
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
  internal_b312c2da0affb0ff73dd81300a784e62::streamcollidesweepthermalizeddoubleprecisionavx_streamcollidesweepthermalizeddoubleprecisionavx(_data_force, _data_pdfs, _data_pdfs_tmp, _size_force_0, _size_force_1, _size_force_2, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_pdfs_tmp_3, block_offset_0, block_offset_1, block_offset_2, kT, omega_bulk, omega_even, omega_odd, omega_shear, seed, time_step);
  pdfs->swapDataPointers(pdfs_tmp);
}

void StreamCollideSweepThermalizedDoublePrecisionAVX::runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks, const CellInterval &globalCellInterval, cell_idx_t ghostLayers, IBlock *block) {
  if (!this->configured_)
    WALBERLA_ABORT("This Sweep contains a configure function that needs to be called manually")

  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto force = block->getData<field::GhostLayerField<double, 3>>(forceID);
  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);
  field::GhostLayerField<double, 19> *pdfs_tmp;
  {
    if (cache_pdfs_.find(block) == cache_pdfs_.end()) {
      pdfs_tmp = pdfs->cloneUninitialized();
      cache_pdfs_[block] = pdfs_tmp;
    } else {
      pdfs_tmp = cache_pdfs_[block];
    }
  }

  auto &seed = this->seed_;
  auto &kT = this->kT_;
  auto &omega_bulk = this->omega_bulk_;
  auto &omega_even = this->omega_even_;
  auto &omega_shear = this->omega_shear_;
  auto &omega_odd = this->omega_odd_;
  auto &block_offset_2 = this->block_offset_2_;
  auto &time_step = this->time_step_;
  auto &block_offset_0 = this->block_offset_0_;
  auto &block_offset_1 = this->block_offset_1_;
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(force->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(force->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(force->nrOfGhostLayers()))
  double *RESTRICT const _data_force = force->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)force->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs->nrOfGhostLayers()))
  double *RESTRICT const _data_pdfs = pdfs->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx)
  WALBERLA_ASSERT_EQUAL((uintptr_t)pdfs->dataAt(0, 0, 0, 0) % 32, 0)
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()))
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_tmp->nrOfGhostLayers()))
  double *RESTRICT _data_pdfs_tmp = pdfs_tmp->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
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
  internal_b312c2da0affb0ff73dd81300a784e62::streamcollidesweepthermalizeddoubleprecisionavx_streamcollidesweepthermalizeddoubleprecisionavx(_data_force, _data_pdfs, _data_pdfs_tmp, _size_force_0, _size_force_1, _size_force_2, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, _stride_pdfs_tmp_1, _stride_pdfs_tmp_2, _stride_pdfs_tmp_3, block_offset_0, block_offset_1, block_offset_2, kT, omega_bulk, omega_even, omega_odd, omega_shear, seed, time_step);
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
