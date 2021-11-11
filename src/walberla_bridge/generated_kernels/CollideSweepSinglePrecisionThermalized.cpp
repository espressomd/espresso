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
//! \\file CollideSweepSinglePrecisionThermalized.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

#include <cmath>

#include "CollideSweepSinglePrecisionThermalized.h"
#include "core/DataTypes.h"
#include "core/Macros.h"

#include "philox_rand.h"

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

namespace internal_collidesweepsingleprecisionthermalized {
static FUNC_PREFIX void collidesweepsingleprecisionthermalized(
    float *RESTRICT const _data_force, float *RESTRICT _data_pdfs,
    int64_t const _size_force_0, int64_t const _size_force_1,
    int64_t const _size_force_2, int64_t const _stride_force_0,
    int64_t const _stride_force_1, int64_t const _stride_force_2,
    int64_t const _stride_force_3, int64_t const _stride_pdfs_0,
    int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
    int64_t const _stride_pdfs_3, uint32_t block_offset_0,
    uint32_t block_offset_1, uint32_t block_offset_2, float kT,
    float omega_bulk, float omega_even, float omega_odd, float omega_shear,
    uint32_t seed, uint32_t time_step) {
  const float xi_25 = -omega_bulk;
  const float xi_36 = -omega_shear;
  const float xi_37 = xi_36 + 2.0f;
  const float xi_38 = xi_37 * 0.5f;
  const float xi_43 = xi_37 * 0.0833333333333333f;
  const float xi_48 = xi_37 * 0.166666666666667f;
  const float xi_58 = xi_37 * 0.25f;
  const float xi_63 = xi_37 * 0.0416666666666667f;
  const float xi_90 = 2.4494897427831779;
  const float xi_115 = omega_odd * 0.25f;
  const float xi_131 = omega_odd * 0.0833333333333333f;
  const float xi_196 = omega_shear * 0.25f;
  const float xi_211 = omega_odd * 0.0416666666666667f;
  const float xi_213 = omega_odd * 0.125f;
  const int64_t rr_0 = 0.0f;
  const float xi_120 = rr_0 * 0.166666666666667f;
  const float xi_186 = rr_0 * 0.0833333333333333f;
  for (int64_t ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1) {
    float *RESTRICT _data_pdfs_20_311 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_31 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_316 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_36 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_31 =
        _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    float *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_34 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_33 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_315 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_32 =
        _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    float *RESTRICT _data_pdfs_20_312 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_35 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_32 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_317 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_318 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    float *RESTRICT _data_pdfs_20_314 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_313 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    float *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1) {
      float *RESTRICT _data_pdfs_20_311_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_311;
      float *RESTRICT _data_pdfs_20_31_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_31;
      float *RESTRICT _data_pdfs_20_316_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_316;
      float *RESTRICT _data_pdfs_20_36_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_36;
      float *RESTRICT _data_force_20_31_10 =
          _stride_force_1 * ctr_1 + _data_force_20_31;
      float *RESTRICT _data_pdfs_20_38_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_38;
      float *RESTRICT _data_pdfs_20_310_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_310;
      float *RESTRICT _data_pdfs_20_34_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      float *RESTRICT _data_pdfs_20_33_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      float *RESTRICT _data_pdfs_20_315_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_315;
      float *RESTRICT _data_force_20_32_10 =
          _stride_force_1 * ctr_1 + _data_force_20_32;
      float *RESTRICT _data_pdfs_20_312_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_312;
      float *RESTRICT _data_pdfs_20_35_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_35;
      float *RESTRICT _data_pdfs_20_32_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_32;
      float *RESTRICT _data_pdfs_20_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_317;
      float *RESTRICT _data_pdfs_20_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_318;
      float *RESTRICT _data_force_20_30_10 =
          _stride_force_1 * ctr_1 + _data_force_20_30;
      float *RESTRICT _data_pdfs_20_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_314;
      float *RESTRICT _data_pdfs_20_39_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_39;
      float *RESTRICT _data_pdfs_20_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_313;
      float *RESTRICT _data_pdfs_20_30_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      float *RESTRICT _data_pdfs_20_37_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_37;
      for (int64_t ctr_0 = 0; ctr_0 < _size_force_0; ctr_0 += 1) {
        const float xi_248 = _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0];
        const float xi_249 = _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0];
        const float xi_250 = _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0];
        const float xi_251 = _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0];
        const float xi_252 = _data_force_20_31_10[_stride_force_0 * ctr_0];
        const float xi_253 = _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0];
        const float xi_254 = _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0];
        const float xi_255 = _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0];
        const float xi_256 = _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0];
        const float xi_257 = _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0];
        const float xi_258 = _data_force_20_32_10[_stride_force_0 * ctr_0];
        const float xi_259 = _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0];
        const float xi_260 = _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0];
        const float xi_261 = _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0];
        const float xi_262 = _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0];
        const float xi_263 = _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0];
        const float xi_264 = _data_force_20_30_10[_stride_force_0 * ctr_0];
        const float xi_265 = _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0];
        const float xi_266 = _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0];
        const float xi_267 = _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0];
        const float xi_268 = _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0];
        const float xi_269 = _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0];

        double random_7_0;
        double random_7_1;
        philox_double2(time_step, block_offset_0 + ctr_0,
                       block_offset_1 + ctr_1, block_offset_2 + ctr_2, 7, seed,
                       random_7_0, random_7_1);

        double random_6_0;
        double random_6_1;
        philox_double2(time_step, block_offset_0 + ctr_0,
                       block_offset_1 + ctr_1, block_offset_2 + ctr_2, 6, seed,
                       random_6_0, random_6_1);

        double random_5_0;
        double random_5_1;
        philox_double2(time_step, block_offset_0 + ctr_0,
                       block_offset_1 + ctr_1, block_offset_2 + ctr_2, 5, seed,
                       random_5_0, random_5_1);

        double random_4_0;
        double random_4_1;
        philox_double2(time_step, block_offset_0 + ctr_0,
                       block_offset_1 + ctr_1, block_offset_2 + ctr_2, 4, seed,
                       random_4_0, random_4_1);

        double random_3_0;
        double random_3_1;
        philox_double2(time_step, block_offset_0 + ctr_0,
                       block_offset_1 + ctr_1, block_offset_2 + ctr_2, 3, seed,
                       random_3_0, random_3_1);

        double random_2_0;
        double random_2_1;
        philox_double2(time_step, block_offset_0 + ctr_0,
                       block_offset_1 + ctr_1, block_offset_2 + ctr_2, 2, seed,
                       random_2_0, random_2_1);

        double random_1_0;
        double random_1_1;
        philox_double2(time_step, block_offset_0 + ctr_0,
                       block_offset_1 + ctr_1, block_offset_2 + ctr_2, 1, seed,
                       random_1_0, random_1_1);

        double random_0_0;
        double random_0_1;
        philox_double2(time_step, block_offset_0 + ctr_0,
                       block_offset_1 + ctr_1, block_offset_2 + ctr_2, 0, seed,
                       random_0_0, random_0_1);

        const float xi_0 = xi_263 + xi_265;
        const float xi_1 = xi_0 + xi_255;
        const float xi_2 = xi_248 + xi_249 + xi_257;
        const float xi_3 = xi_259 + xi_260;
        const float xi_4 = xi_256 + xi_266;
        const float xi_5 = xi_250 + xi_261;
        const float xi_6 = xi_251 + xi_262;
        const float xi_8 = -xi_266;
        const float xi_9 = -xi_269 + xi_8;
        const float xi_10 = -xi_262;
        const float xi_11 = -xi_267;
        const float xi_12 = -xi_256;
        const float xi_13 = xi_10 + xi_11 + xi_12;
        const float xi_14 = -xi_261;
        const float xi_15 = -xi_254;
        const float xi_16 = xi_14 + xi_15;
        const float xi_17 = -xi_250;
        const float xi_18 = -xi_259;
        const float xi_19 = xi_17 + xi_18;
        const float xi_20 = -xi_263;
        const float xi_21 = xi_10 + xi_20;
        const float xi_22 = -xi_257;
        const float xi_23 = -xi_251;
        const float xi_24 = xi_17 + xi_22 + xi_23 + xi_248;
        const float xi_42 = xi_252 * 0.166666666666667f;
        const float xi_50 = xi_264 * 0.166666666666667f;
        const float xi_54 = xi_258 * 0.166666666666667f;
        const float xi_57 = xi_252 * 0.5f;
        const float xi_61 = xi_264 * 0.0833333333333333f;
        const float xi_65 = xi_252 * 0.0833333333333333f;
        const float xi_75 = xi_258 * 0.0833333333333333f;
        const float xi_93 = -xi_268;
        const float xi_94 = xi_251 * 3.0f + xi_260 * 3.0f + xi_93;
        const float xi_95 =
            omega_even *
            (xi_248 * -3.0f + xi_249 * 3.0f + xi_250 * -3.0f + xi_257 * -3.0f +
             xi_259 * -3.0f + xi_261 * 3.0f + xi_94);
        const float xi_96 =
            xi_248 * 2.0f + xi_250 * 2.0f + xi_257 * 2.0f + xi_259 * 2.0f;
        const float xi_97 = xi_255 * 5.0f + xi_256 * 5.0f + xi_96;
        const float xi_98 =
            omega_even *
            (xi_249 * -2.0f + xi_261 * -2.0f + xi_262 * -5.0f + xi_263 * -5.0f +
             xi_265 * -5.0f + xi_267 * -5.0f + xi_94 + xi_97);
        const float xi_101 = -xi_248;
        const float xi_102 = xi_101 + xi_18;
        const float xi_103 = -xi_253;
        const float xi_106 = -xi_265;
        const float xi_107 = xi_106 + xi_11 + xi_15 + xi_21;
        const float xi_109 = xi_267 * 2.0f;
        const float xi_110 = xi_265 * 2.0f;
        const float xi_111 = xi_262 * 2.0f + xi_263 * 2.0f;
        const float xi_112 =
            omega_even *
            (xi_109 + xi_110 + xi_111 + xi_249 * 5.0f + xi_251 * -4.0f +
             xi_253 * -7.0f + xi_254 * -7.0f + xi_260 * -4.0f + xi_261 * 5.0f +
             xi_266 * -7.0f + xi_269 * -7.0f + xi_93 + xi_97);
        const float xi_113 = xi_101 + xi_259;
        const float xi_114 = xi_113 + xi_14 + xi_22 + xi_249 + xi_250;
        const float xi_116 = xi_114 * xi_115;
        const float xi_118 = xi_103 + xi_254;
        const double xi_122 = random_5_1 - 0.5f;
        const float xi_127 = xi_269 * 2.0f;
        const float xi_128 = xi_254 * 2.0f;
        const float xi_129 = xi_253 * -2.0f + xi_266 * 2.0f;
        const float xi_130 = -xi_127 + xi_128 + xi_129 + xi_14 + xi_19 + xi_2;
        const float xi_132 = xi_130 * xi_131;
        const double xi_133 = random_3_0 - 0.5f;
        const double xi_138 = random_0_1 - 0.5f;
        const float xi_142 = xi_262 + xi_267;
        const float xi_156 = xi_106 + xi_267;
        const float xi_157 = xi_12 + xi_156 + xi_20 + xi_255 + xi_262;
        const float xi_158 = xi_115 * xi_157;
        const double xi_159 = random_4_1 - 0.5f;
        const float xi_161 = xi_1 + xi_127 - xi_128 + xi_129 + xi_13;
        const float xi_162 = xi_131 * xi_161;
        const double xi_163 = random_4_0 - 0.5f;
        const float xi_168 = xi_250 + xi_257;
        const float xi_169 = xi_102 + xi_168 + xi_23 + xi_260;
        const float xi_170 = xi_115 * xi_169;
        const double xi_173 = random_5_0 - 0.5f;
        const float xi_175 = -xi_109 - xi_110 + xi_111 + xi_24 + xi_3;
        const float xi_176 = xi_131 * xi_175;
        const double xi_177 = random_3_1 - 0.5f;
        const float xi_184 = xi_112 * 0.0138888888888889f;
        const float xi_205 = xi_98 * -0.00714285714285714f;
        const float xi_207 = xi_95 * 0.025f;
        const float xi_212 = xi_175 * xi_211;
        const float xi_214 = xi_169 * xi_213;
        const float xi_223 = xi_130 * xi_211;
        const float xi_224 = xi_114 * xi_213;
        const float xi_232 = xi_98 * 0.0178571428571429f;
        const float xi_238 = xi_157 * xi_213;
        const float xi_239 = xi_161 * xi_211;
        const float vel0Term = xi_1 + xi_253 + xi_254;
        const float vel1Term = xi_2 + xi_269;
        const float vel2Term = xi_267 + xi_3;
        const float rho =
            vel0Term + vel1Term + vel2Term + xi_268 + xi_4 + xi_5 + xi_6;
        const float xi_7 = 1 / (rho);
        const float xi_86 = kT * rho;
        const float xi_87 = sqrt(
            xi_86 * (-((-omega_even + 1.0f) * (-omega_even + 1.0f)) + 1.0f));
        const float xi_88 = xi_87 * (random_6_0 - 0.5f) * 3.7416573867739413;
        const float xi_89 = xi_87 * (random_7_0 - 0.5f) * 5.4772255750516612;
        const float xi_91 =
            xi_90 * sqrt(xi_86 * (-((xi_25 + 1.0f) * (xi_25 + 1.0f)) + 1.0f)) *
            (random_2_1 - 0.5f);
        const float xi_92 = xi_87 * (random_6_1 - 0.5f) * 8.3666002653407556;
        const float xi_123 =
            sqrt(xi_86 * (-((-omega_odd + 1.0f) * (-omega_odd + 1.0f)) + 1.0f));
        const float xi_124 = xi_123 * 1.4142135623730951;
        const float xi_125 = xi_124 * 0.5f;
        const float xi_126 = xi_122 * xi_125;
        const float xi_134 = xi_123 * xi_90;
        const float xi_135 = xi_134 * 0.166666666666667f;
        const float xi_136 = xi_133 * xi_135;
        const float xi_137 = -xi_132 - xi_136;
        const float xi_139 =
            sqrt(xi_86 * (-((xi_36 + 1.0f) * (xi_36 + 1.0f)) + 1.0f));
        const float xi_140 = xi_139 * 0.5f;
        const float xi_141 = xi_138 * xi_140;
        const float xi_146 =
            xi_112 * -0.0198412698412698f + xi_88 * -0.119047619047619f;
        const float xi_148 = xi_139 * (random_0_0 - 0.5f) * 1.7320508075688772;
        const float xi_152 = xi_132 + xi_136;
        const float xi_160 = xi_125 * xi_159;
        const float xi_164 = xi_135 * xi_163;
        const float xi_165 = xi_162 + xi_164;
        const float xi_167 = -xi_162 - xi_164;
        const float xi_174 = xi_125 * xi_173;
        const float xi_178 = xi_135 * xi_177;
        const float xi_179 = -xi_176 - xi_178;
        const float xi_181 = xi_176 + xi_178;
        const float xi_182 = xi_138 * xi_139 * 0.25f;
        const float xi_185 = xi_88 * 0.0833333333333333f;
        const float xi_195 = xi_140 * (random_1_0 - 0.5f);
        const float xi_204 = xi_140 * (random_2_0 - 0.5f);
        const float xi_208 = xi_92 * -0.0142857142857143f;
        const float xi_209 = xi_89 * 0.05f;
        const float xi_215 = xi_134 * 0.0833333333333333f;
        const float xi_216 = xi_177 * xi_215;
        const float xi_217 = xi_124 * 0.25f;
        const float xi_218 = xi_173 * xi_217;
        const float xi_220 =
            xi_112 * -0.00396825396825397f + xi_88 * -0.0238095238095238f;
        const float xi_225 = xi_133 * xi_215;
        const float xi_226 = xi_122 * xi_217;
        const float xi_230 = -xi_182;
        const float xi_233 = xi_92 * 0.0357142857142857f;
        const float xi_235 = xi_140 * (random_1_1 - 0.5f);
        const float xi_240 = xi_159 * xi_217;
        const float xi_241 = xi_163 * xi_215;
        const float u_0 = xi_7 * (vel0Term + xi_13 + xi_9);
        const float xi_26 = u_0 * xi_264;
        const float xi_27 = xi_26 * 0.333333333333333f;
        const float xi_33 = -xi_27;
        const float xi_99 = rho * (u_0 * u_0);
        const float xi_153 = rho * u_0;
        const float xi_154 = -vel0Term + xi_142 + xi_153 + xi_269 + xi_4;
        const float xi_155 = xi_120 * xi_154;
        const float xi_191 = xi_154 * xi_186;
        const float u_1 = xi_7 * (vel1Term + xi_16 + xi_19 + xi_253 + xi_8);
        const float xi_28 = u_1 * xi_252;
        const float xi_29 = xi_28 * 0.333333333333333f;
        const float xi_34 = -xi_29;
        const float xi_56 = u_1 * 0.5f;
        const float xi_59 = xi_58 * (u_0 * xi_57 + xi_264 * xi_56);
        const float xi_60 = -xi_59;
        const float xi_104 = rho * (u_1 * u_1);
        const float xi_105 = xi_103 + xi_104 + xi_9;
        const float xi_117 = rho * u_1;
        const float xi_119 =
            -vel1Term + xi_117 + xi_118 + xi_259 + xi_266 + xi_5;
        const float xi_121 = xi_119 * xi_120;
        const float xi_187 = xi_119 * xi_186;
        const float xi_197 = xi_196 * (u_0 * xi_117 + xi_118 + xi_269 + xi_8);
        const float xi_198 = -xi_195 - xi_197;
        const float xi_199 = xi_195 + xi_197;
        const float u_2 = xi_7 * (vel2Term + xi_21 + xi_24 + xi_265);
        const float xi_30 = u_2 * xi_258;
        const float xi_31 = xi_30 * 0.333333333333333f;
        const float xi_32 = (xi_25 + 2.0f) * (xi_27 + xi_29 + xi_31);
        const float xi_35 = xi_30 * 0.666666666666667f + xi_33 + xi_34;
        const float xi_39 = -xi_31;
        const float xi_40 = xi_28 * 0.666666666666667f + xi_33 + xi_39;
        const float xi_41 = xi_26 * 0.666666666666667f + xi_34 + xi_39;
        const float xi_44 = xi_35 * xi_43;
        const float xi_45 = -xi_44;
        const float xi_46 = xi_41 * xi_43;
        const float xi_47 = -xi_46;
        const float xi_49 = xi_40 * xi_48 + xi_45 + xi_47;
        const float xi_51 = xi_40 * xi_43;
        const float xi_52 = -xi_51;
        const float xi_53 = xi_41 * xi_48 + xi_45 + xi_52;
        const float xi_55 = xi_35 * xi_48 + xi_47 + xi_52;
        const float xi_62 = xi_46 - xi_61;
        const float xi_64 = -xi_35 * xi_63;
        const float xi_66 = xi_32 * 0.125f;
        const float xi_67 = xi_51 + xi_66;
        const float xi_68 = xi_65 + xi_67;
        const float xi_69 = xi_64 + xi_68;
        const float xi_70 = xi_46 + xi_61;
        const float xi_71 = -xi_65 + xi_67;
        const float xi_72 = xi_64 + xi_71;
        const float xi_73 = xi_58 * (u_2 * xi_57 + xi_258 * xi_56);
        const float xi_74 = -xi_41 * xi_63;
        const float xi_76 = xi_44 + xi_75;
        const float xi_77 = xi_74 + xi_76;
        const float xi_78 = -xi_73;
        const float xi_79 = xi_58 * (u_0 * xi_258 * 0.5f + u_2 * xi_264 * 0.5f);
        const float xi_80 = -xi_79;
        const float xi_81 = -xi_40 * xi_63;
        const float xi_82 = xi_66 + xi_76 + xi_81;
        const float xi_83 = xi_44 - xi_75;
        const float xi_84 = xi_74 + xi_83;
        const float xi_85 = xi_66 + xi_81 + xi_83;
        const float xi_100 = rho * (u_2 * u_2);
        const float xi_108 = omega_bulk * (xi_100 + xi_102 + xi_105 + xi_107 +
                                           xi_17 + xi_22 + xi_268 + xi_99);
        const float xi_143 = -xi_100 + xi_251 + xi_260;
        const float xi_144 =
            omega_shear * (xi_0 + xi_105 + xi_142 + xi_143 + xi_16 - xi_249);
        const float xi_145 = xi_144 * 0.125f;
        const float xi_147 =
            omega_shear *
            (xi_103 - xi_104 + xi_107 + xi_143 + xi_249 + xi_255 * -2.0f +
             xi_256 * -2.0f + xi_261 + xi_9 + xi_96 + xi_99 * 2.0f);
        const float xi_149 =
            xi_147 * -0.0416666666666667f + xi_148 * -0.166666666666667f;
        const float xi_150 = xi_149 + xi_89 * -0.1f + xi_95 * -0.05f;
        const float xi_151 = xi_141 + xi_145 + xi_146 + xi_150 +
                             xi_92 * 0.0285714285714286f +
                             xi_98 * 0.0142857142857143f;
        const float xi_166 = xi_146 + xi_147 * 0.0833333333333333f +
                             xi_148 * 0.333333333333333f +
                             xi_92 * -0.0714285714285714f +
                             xi_98 * -0.0357142857142857f;
        const float xi_171 =
            rho * u_2 - vel2Term + xi_101 + xi_106 + xi_168 + xi_263 + xi_6;
        const float xi_172 = xi_120 * xi_171;
        const float xi_180 = xi_112 * 0.0158730158730159f - xi_141 - xi_145 +
                             xi_150 + xi_88 * 0.0952380952380952f +
                             xi_92 * -0.0428571428571429f +
                             xi_98 * -0.0214285714285714f;
        const float xi_183 = xi_144 * 0.0625f;
        const float xi_188 =
            xi_108 * 0.0416666666666667f + xi_91 * 0.0833333333333333f;
        const float xi_189 = xi_187 + xi_188;
        const float xi_190 =
            xi_152 + xi_182 + xi_183 + xi_184 + xi_185 + xi_189;
        const float xi_192 =
            xi_147 * 0.0208333333333333f + xi_148 * 0.0833333333333333f;
        const float xi_193 = -xi_191 + xi_192;
        const float xi_194 = xi_167 + xi_193;
        const float xi_200 = xi_191 + xi_192;
        const float xi_201 = xi_165 + xi_200;
        const float xi_202 = -xi_187 + xi_188;
        const float xi_203 =
            xi_137 + xi_182 + xi_183 + xi_184 + xi_185 + xi_202;
        const float xi_206 = xi_196 * (u_2 * xi_117 + xi_113 + xi_17 + xi_257);
        const float xi_210 =
            xi_149 + xi_204 + xi_205 + xi_206 + xi_207 + xi_208 + xi_209;
        const float xi_219 = xi_171 * xi_186;
        const float xi_221 = xi_219 + xi_220;
        const float xi_222 = -xi_212 + xi_214 - xi_216 + xi_218 + xi_221;
        const float xi_227 = xi_189 - xi_223 + xi_224 - xi_225 + xi_226;
        const float xi_228 = xi_202 + xi_223 - xi_224 + xi_225 - xi_226;
        const float xi_229 =
            xi_149 - xi_204 + xi_205 - xi_206 + xi_207 + xi_208 + xi_209;
        const float xi_231 = -xi_183;
        const float xi_234 =
            xi_181 + xi_188 + xi_221 + xi_230 + xi_231 + xi_232 + xi_233;
        const float xi_236 = xi_196 * (u_2 * xi_153 + xi_10 + xi_156 + xi_263);
        const float xi_237 = -xi_235 - xi_236;
        const float xi_242 = xi_193 - xi_238 + xi_239 - xi_240 + xi_241;
        const float xi_243 = xi_235 + xi_236;
        const float xi_244 = xi_200 + xi_238 - xi_239 + xi_240 - xi_241;
        const float xi_245 = -xi_219 + xi_220;
        const float xi_246 = xi_212 - xi_214 + xi_216 - xi_218 + xi_245;
        const float xi_247 =
            xi_179 + xi_188 + xi_230 + xi_231 + xi_232 + xi_233 + xi_245;
        const float forceTerm_0 =
            xi_32 * -1.5f - xi_35 * xi_38 - xi_38 * xi_40 - xi_38 * xi_41;
        const float forceTerm_1 = xi_42 + xi_49;
        const float forceTerm_2 = -xi_42 + xi_49;
        const float forceTerm_3 = -xi_50 + xi_53;
        const float forceTerm_4 = xi_50 + xi_53;
        const float forceTerm_5 = xi_54 + xi_55;
        const float forceTerm_6 = -xi_54 + xi_55;
        const float forceTerm_7 = xi_60 + xi_62 + xi_69;
        const float forceTerm_8 = xi_59 + xi_69 + xi_70;
        const float forceTerm_9 = xi_59 + xi_62 + xi_72;
        const float forceTerm_10 = xi_60 + xi_70 + xi_72;
        const float forceTerm_11 = xi_68 + xi_73 + xi_77;
        const float forceTerm_12 = xi_71 + xi_77 + xi_78;
        const float forceTerm_13 = xi_62 + xi_80 + xi_82;
        const float forceTerm_14 = xi_70 + xi_79 + xi_82;
        const float forceTerm_15 = xi_68 + xi_78 + xi_84;
        const float forceTerm_16 = xi_71 + xi_73 + xi_84;
        const float forceTerm_17 = xi_62 + xi_79 + xi_85;
        const float forceTerm_18 = xi_70 + xi_80 + xi_85;
        _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_0 + xi_108 * -0.5f + xi_112 * 0.0238095238095238f +
            xi_268 + xi_88 * 0.142857142857143f + xi_89 * 0.2f - xi_91 +
            xi_92 * 0.0857142857142857f + xi_95 * 0.1f +
            xi_98 * 0.0428571428571429f;
        _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_1 - xi_116 + xi_121 - xi_126 + xi_137 + xi_151 + xi_249;
        _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_2 + xi_116 - xi_121 + xi_126 + xi_151 + xi_152 + xi_261;
        _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_3 - xi_155 + xi_158 + xi_160 + xi_165 + xi_166 + xi_256;
        _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_4 + xi_155 - xi_158 - xi_160 + xi_166 + xi_167 + xi_255;
        _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_5 - xi_170 + xi_172 - xi_174 + xi_179 + xi_180 + xi_260;
        _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_6 + xi_170 - xi_172 + xi_174 + xi_180 + xi_181 + xi_251;
        _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_7 + xi_190 + xi_194 + xi_198 + xi_269;
        _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_8 + xi_190 + xi_199 + xi_201 + xi_253;
        _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_9 + xi_194 + xi_199 + xi_203 + xi_266;
        _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_10 + xi_198 + xi_201 + xi_203 + xi_254;
        _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_11 + xi_210 + xi_222 + xi_227 + xi_248;
        _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_12 + xi_222 + xi_228 + xi_229 + xi_259;
        _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_13 + xi_234 + xi_237 + xi_242 + xi_267;
        _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_14 + xi_234 + xi_243 + xi_244 + xi_265;
        _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_15 + xi_227 + xi_229 + xi_246 + xi_257;
        _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_16 + xi_210 + xi_228 + xi_246 + xi_250;
        _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_17 + xi_242 + xi_243 + xi_247 + xi_262;
        _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_18 + xi_237 + xi_244 + xi_247 + xi_263;
      }
    }
  }
}
} // namespace internal_collidesweepsingleprecisionthermalized

void CollideSweepSinglePrecisionThermalized::operator()(IBlock *block) {
  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);
  auto force = block->getData<field::GhostLayerField<float, 3>>(forceID);

  auto &omega_odd = this->omega_odd_;
  auto &omega_shear = this->omega_shear_;
  auto &seed = this->seed_;
  auto &kT = this->kT_;
  auto block_offset_2 = this->block_offset_2_;
  auto &time_step = this->time_step_;
  auto block_offset_1 = this->block_offset_1_;
  auto &omega_even = this->omega_even_;
  auto &omega_bulk = this->omega_bulk_;
  auto block_offset_0 = this->block_offset_0_;
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
  const int64_t _stride_force_0 = int64_t(force->xStride());
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  internal_collidesweepsingleprecisionthermalized::
      collidesweepsingleprecisionthermalized(
          _data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2,
          _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3,
          _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3,
          block_offset_0, block_offset_1, block_offset_2, kT, omega_bulk,
          omega_even, omega_odd, omega_shear, seed, time_step);
}

void CollideSweepSinglePrecisionThermalized::runOnCellInterval(
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

  auto &omega_odd = this->omega_odd_;
  auto &omega_shear = this->omega_shear_;
  auto &seed = this->seed_;
  auto &kT = this->kT_;
  auto block_offset_2 = this->block_offset_2_;
  auto &time_step = this->time_step_;
  auto block_offset_1 = this->block_offset_1_;
  auto &omega_even = this->omega_even_;
  auto &omega_bulk = this->omega_bulk_;
  auto block_offset_0 = this->block_offset_0_;
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
  const int64_t _stride_force_0 = int64_t(force->xStride());
  const int64_t _stride_force_1 = int64_t(force->yStride());
  const int64_t _stride_force_2 = int64_t(force->zStride());
  const int64_t _stride_force_3 = int64_t(1 * int64_t(force->fStride()));
  const int64_t _stride_pdfs_0 = int64_t(pdfs->xStride());
  const int64_t _stride_pdfs_1 = int64_t(pdfs->yStride());
  const int64_t _stride_pdfs_2 = int64_t(pdfs->zStride());
  const int64_t _stride_pdfs_3 = int64_t(1 * int64_t(pdfs->fStride()));
  internal_collidesweepsingleprecisionthermalized::
      collidesweepsingleprecisionthermalized(
          _data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2,
          _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3,
          _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3,
          block_offset_0, block_offset_1, block_offset_2, kT, omega_bulk,
          omega_even, omega_odd, omega_shear, seed, time_step);
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
