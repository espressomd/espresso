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
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning push
#pragma warning(disable : 1599)
#endif

using namespace std;

namespace walberla {
namespace pystencils {

namespace internal_69764eed2d0964e29e3b97d1054b4693 {
static FUNC_PREFIX void
collidesweepsingleprecisionthermalized_collidesweepsingleprecisionthermalized(
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
  const float xi_29 = omega_bulk * 0.50000000000000000f;
  const float xi_56 = omega_shear * 0.041666666666666667f;
  const float xi_61 = omega_bulk * 0.041666666666666667f;
  const float xi_72 = omega_shear * 0.12500000000000000f;
  const float xi_110 = 2.4494897427831781f;
  const float xi_135 = omega_odd * 0.25000000000000000f;
  const float xi_151 = omega_odd * 0.083333333333333333f;
  const float xi_216 = omega_shear * 0.25000000000000000f;
  const float xi_231 = omega_odd * 0.041666666666666667f;
  const float xi_233 = omega_odd * 0.12500000000000000f;
  const float rr_0 = 0.0f;
  const float xi_54 = rr_0 * 0.041666666666666667f;
  const float xi_140 = rr_0 * 0.16666666666666667f;
  const float xi_206 = rr_0 * 0.083333333333333333f;
  for (int64_t ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1) {
    float *RESTRICT _data_pdfs_20_31 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_34 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    float *RESTRICT _data_pdfs_20_316 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_314 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_311 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_317 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    float *RESTRICT _data_pdfs_20_318 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_33 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_315 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_312 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_31 =
        _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    float *RESTRICT _data_force_20_32 =
        _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    float *RESTRICT _data_pdfs_20_32 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_35 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_313 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_36 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1) {
      float *RESTRICT _data_pdfs_20_31_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_31;
      float *RESTRICT _data_pdfs_20_39_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_39;
      float *RESTRICT _data_pdfs_20_34_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      float *RESTRICT _data_force_20_30_10 =
          _stride_force_1 * ctr_1 + _data_force_20_30;
      float *RESTRICT _data_pdfs_20_316_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_316;
      float *RESTRICT _data_pdfs_20_37_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_37;
      float *RESTRICT _data_pdfs_20_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_314;
      float *RESTRICT _data_pdfs_20_311_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_311;
      float *RESTRICT _data_pdfs_20_38_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_38;
      float *RESTRICT _data_pdfs_20_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_317;
      float *RESTRICT _data_pdfs_20_30_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      float *RESTRICT _data_pdfs_20_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_318;
      float *RESTRICT _data_pdfs_20_33_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      float *RESTRICT _data_pdfs_20_315_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_315;
      float *RESTRICT _data_pdfs_20_312_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_312;
      float *RESTRICT _data_force_20_31_10 =
          _stride_force_1 * ctr_1 + _data_force_20_31;
      float *RESTRICT _data_force_20_32_10 =
          _stride_force_1 * ctr_1 + _data_force_20_32;
      float *RESTRICT _data_pdfs_20_32_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_32;
      float *RESTRICT _data_pdfs_20_35_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_35;
      float *RESTRICT _data_pdfs_20_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_313;
      float *RESTRICT _data_pdfs_20_310_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_310;
      float *RESTRICT _data_pdfs_20_36_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_36;
      for (int64_t ctr_0 = 0; ctr_0 < _size_force_0; ctr_0 += 1) {
        const float xi_268 = _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0];
        const float xi_269 = _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0];
        const float xi_270 = _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0];
        const float xi_271 = _data_force_20_30_10[_stride_force_0 * ctr_0];
        const float xi_272 = _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0];
        const float xi_273 = _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0];
        const float xi_274 = _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0];
        const float xi_275 = _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0];
        const float xi_276 = _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0];
        const float xi_277 = _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0];
        const float xi_278 = _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0];
        const float xi_279 = _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0];
        const float xi_280 = _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0];
        const float xi_281 = _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0];
        const float xi_282 = _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0];
        const float xi_283 = _data_force_20_31_10[_stride_force_0 * ctr_0];
        const float xi_284 = _data_force_20_32_10[_stride_force_0 * ctr_0];
        const float xi_285 = _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0];
        const float xi_286 = _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0];
        const float xi_287 = _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0];
        const float xi_288 = _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0];
        const float xi_289 = _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0];

        float random_3_0;
        float random_3_1;
        float random_3_2;
        float random_3_3;
        philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1,
                      block_offset_2 + ctr_2, 3, seed, random_3_0, random_3_1,
                      random_3_2, random_3_3);

        float random_2_0;
        float random_2_1;
        float random_2_2;
        float random_2_3;
        philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1,
                      block_offset_2 + ctr_2, 2, seed, random_2_0, random_2_1,
                      random_2_2, random_2_3);

        float random_1_0;
        float random_1_1;
        float random_1_2;
        float random_1_3;
        philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1,
                      block_offset_2 + ctr_2, 1, seed, random_1_0, random_1_1,
                      random_1_2, random_1_3);

        float random_0_0;
        float random_0_1;
        float random_0_2;
        float random_0_3;
        philox_float4(time_step, block_offset_0 + ctr_0, block_offset_1 + ctr_1,
                      block_offset_2 + ctr_2, 0, seed, random_0_0, random_0_1,
                      random_0_2, random_0_3);

        const float xi_0 = xi_274 + xi_279;
        const float xi_1 = xi_0 + xi_270;
        const float xi_2 = xi_268 + xi_275 + xi_281;
        const float xi_3 = xi_282 + xi_286;
        const float xi_4 = xi_269 + xi_280;
        const float xi_5 = xi_272 + xi_285;
        const float xi_6 = xi_277 + xi_289;
        const float xi_9 = -xi_269;
        const float xi_10 = -xi_273 + xi_9;
        const float xi_11 = -xi_277;
        const float xi_12 = -xi_287;
        const float xi_13 = -xi_280;
        const float xi_14 = xi_11 + xi_12 + xi_13;
        const float xi_15 = -xi_285;
        const float xi_16 = -xi_288;
        const float xi_17 = xi_15 + xi_16;
        const float xi_18 = -xi_272;
        const float xi_19 = -xi_282;
        const float xi_20 = xi_18 + xi_19;
        const float xi_21 = -xi_279;
        const float xi_22 = xi_11 + xi_21;
        const float xi_23 = -xi_281;
        const float xi_24 = -xi_289;
        const float xi_25 = xi_18 + xi_23 + xi_24 + xi_275;
        const float xi_30 = xi_283 * 0.16666666666666667f;
        const float xi_31 = xi_283 * 0.083333333333333333f;
        const float xi_43 = xi_271 * 0.16666666666666667f;
        const float xi_44 = xi_271 * 0.083333333333333333f;
        const float xi_50 = xi_284 * 0.16666666666666667f;
        const float xi_51 = xi_284 * 0.083333333333333333f;
        const float xi_68 = xi_283 * 0.25000000000000000f;
        const float xi_73 = xi_283 * xi_72;
        const float xi_113 = -xi_278;
        const float xi_114 = xi_113 + xi_286 * 3.0f + xi_289 * 3.0f;
        const float xi_115 =
            omega_even *
            (xi_114 + xi_268 * 3.0f + xi_272 * -3.0f + xi_275 * -3.0f +
             xi_281 * -3.0f + xi_282 * -3.0f + xi_285 * 3.0f);
        const float xi_116 =
            xi_272 * 2.0f + xi_275 * 2.0f + xi_281 * 2.0f + xi_282 * 2.0f;
        const float xi_117 = xi_116 + xi_270 * 5.0f + xi_280 * 5.0f;
        const float xi_118 =
            omega_even *
            (xi_114 + xi_117 + xi_268 * -2.0f + xi_274 * -5.0f +
             xi_277 * -5.0f + xi_279 * -5.0f + xi_285 * -2.0f + xi_287 * -5.0f);
        const float xi_121 = -xi_275;
        const float xi_122 = xi_121 + xi_19;
        const float xi_123 = -xi_276;
        const float xi_126 = -xi_274;
        const float xi_127 = xi_12 + xi_126 + xi_16 + xi_22;
        const float xi_129 = xi_287 * 2.0f;
        const float xi_130 = xi_274 * 2.0f;
        const float xi_131 = xi_277 * 2.0f + xi_279 * 2.0f;
        const float xi_132 =
            omega_even *
            (xi_113 + xi_117 + xi_129 + xi_130 + xi_131 + xi_268 * 5.0f +
             xi_269 * -7.0f + xi_273 * -7.0f + xi_276 * -7.0f + xi_285 * 5.0f +
             xi_286 * -4.0f + xi_288 * -7.0f + xi_289 * -4.0f);
        const float xi_133 = xi_121 + xi_282;
        const float xi_134 = xi_133 + xi_15 + xi_23 + xi_268 + xi_272;
        const float xi_136 = xi_134 * xi_135;
        const float xi_138 = xi_123 + xi_288;
        const float xi_142 = random_2_3 - 0.5f;
        const float xi_147 = xi_273 * 2.0f;
        const float xi_148 = xi_288 * 2.0f;
        const float xi_149 = xi_269 * 2.0f + xi_276 * -2.0f;
        const float xi_150 = -xi_147 + xi_148 + xi_149 + xi_15 + xi_2 + xi_20;
        const float xi_152 = xi_150 * xi_151;
        const float xi_153 = random_1_2 - 0.5f;
        const float xi_158 = random_0_1 - 0.5f;
        const float xi_162 = xi_277 + xi_287;
        const float xi_176 = xi_126 + xi_287;
        const float xi_177 = xi_13 + xi_176 + xi_21 + xi_270 + xi_277;
        const float xi_178 = xi_135 * xi_177;
        const float xi_179 = random_2_1 - 0.5f;
        const float xi_181 = xi_1 + xi_14 + xi_147 - xi_148 + xi_149;
        const float xi_182 = xi_151 * xi_181;
        const float xi_183 = random_2_0 - 0.5f;
        const float xi_188 = xi_272 + xi_281;
        const float xi_189 = xi_122 + xi_188 + xi_24 + xi_286;
        const float xi_190 = xi_135 * xi_189;
        const float xi_193 = random_2_2 - 0.5f;
        const float xi_195 = -xi_129 - xi_130 + xi_131 + xi_25 + xi_3;
        const float xi_196 = xi_151 * xi_195;
        const float xi_197 = random_1_3 - 0.5f;
        const float xi_204 = xi_132 * 0.013888888888888889f;
        const float xi_225 = xi_118 * -0.0071428571428571429f;
        const float xi_227 = xi_115 * 0.025000000000000000f;
        const float xi_232 = xi_195 * xi_231;
        const float xi_234 = xi_189 * xi_233;
        const float xi_243 = xi_150 * xi_231;
        const float xi_244 = xi_134 * xi_233;
        const float xi_252 = xi_118 * 0.017857142857142857f;
        const float xi_258 = xi_177 * xi_233;
        const float xi_259 = xi_181 * xi_231;
        const float xi_32 = rr_0 * xi_31;
        const float xi_45 = rr_0 * xi_44;
        const float xi_52 = rr_0 * xi_51;
        const float xi_55 = xi_271 * xi_54;
        const float xi_60 = xi_283 * xi_54;
        const float xi_82 = xi_284 * xi_54;
        const float vel0Term = xi_1 + xi_276 + xi_288;
        const float vel1Term = xi_2 + xi_273;
        const float vel2Term = xi_287 + xi_3;
        const float rho =
            vel0Term + vel1Term + vel2Term + xi_278 + xi_4 + xi_5 + xi_6;
        const float xi_7 = 1 / (rho);
        const float xi_8 = xi_7 * 0.50000000000000000f;
        const float xi_106 = kT * rho;
        const float xi_107 = sqrt(
            xi_106 * (-((-omega_even + 1.0f) * (-omega_even + 1.0f)) + 1.0f));
        const float xi_108 = xi_107 * (random_3_0 - 0.5f) * 3.7416573867739414f;
        const float xi_109 = xi_107 * (random_3_2 - 0.5f) * 5.4772255750516611f;
        const float xi_111 =
            xi_110 *
            sqrt(xi_106 *
                 (-((-omega_bulk + 1.0f) * (-omega_bulk + 1.0f)) + 1.0f)) *
            (random_1_1 - 0.5f);
        const float xi_112 = xi_107 * (random_3_1 - 0.5f) * 8.3666002653407555f;
        const float xi_143 = sqrt(
            xi_106 * (-((-omega_odd + 1.0f) * (-omega_odd + 1.0f)) + 1.0f));
        const float xi_144 = xi_143 * 1.414213562373095f;
        const float xi_145 = xi_144 * 0.50000000000000000f;
        const float xi_146 = xi_142 * xi_145;
        const float xi_154 = xi_110 * xi_143;
        const float xi_155 = xi_154 * 0.16666666666666667f;
        const float xi_156 = xi_153 * xi_155;
        const float xi_157 = -xi_152 - xi_156;
        const float xi_159 = sqrt(
            xi_106 * (-((-omega_shear + 1.0f) * (-omega_shear + 1.0f)) + 1.0f));
        const float xi_160 = xi_159 * 0.50000000000000000f;
        const float xi_161 = xi_158 * xi_160;
        const float xi_166 =
            xi_108 * -0.11904761904761905f + xi_132 * -0.019841269841269841f;
        const float xi_168 = xi_159 * (random_0_0 - 0.5f) * 1.7320508075688773f;
        const float xi_172 = xi_152 + xi_156;
        const float xi_180 = xi_145 * xi_179;
        const float xi_184 = xi_155 * xi_183;
        const float xi_185 = xi_182 + xi_184;
        const float xi_187 = -xi_182 - xi_184;
        const float xi_194 = xi_145 * xi_193;
        const float xi_198 = xi_155 * xi_197;
        const float xi_199 = -xi_196 - xi_198;
        const float xi_201 = xi_196 + xi_198;
        const float xi_202 = xi_158 * xi_159 * 0.25000000000000000f;
        const float xi_205 = xi_108 * 0.083333333333333333f;
        const float xi_215 = xi_160 * (random_0_2 - 0.5f);
        const float xi_224 = xi_160 * (random_1_0 - 0.5f);
        const float xi_228 = xi_112 * -0.014285714285714286f;
        const float xi_229 = xi_109 * 0.050000000000000000f;
        const float xi_235 = xi_154 * 0.083333333333333333f;
        const float xi_236 = xi_197 * xi_235;
        const float xi_237 = xi_144 * 0.25000000000000000f;
        const float xi_238 = xi_193 * xi_237;
        const float xi_240 =
            xi_108 * -0.023809523809523810f + xi_132 * -0.0039682539682539683f;
        const float xi_245 = xi_153 * xi_235;
        const float xi_246 = xi_142 * xi_237;
        const float xi_250 = -xi_202;
        const float xi_253 = xi_112 * 0.035714285714285714f;
        const float xi_255 = xi_160 * (random_0_3 - 0.5f);
        const float xi_260 = xi_179 * xi_237;
        const float xi_261 = xi_183 * xi_235;
        const float u_0 = xi_271 * xi_8 + xi_7 * (vel0Term + xi_10 + xi_14);
        const float xi_26 = u_0 * xi_271;
        const float xi_38 = xi_26 * 0.16666666666666667f;
        const float xi_39 = xi_26 * 0.083333333333333333f;
        const float xi_40 = omega_shear * xi_39;
        const float xi_41 = -xi_38 + xi_40;
        const float xi_57 = -xi_26 * xi_56 + xi_38;
        const float xi_58 = -xi_44 + xi_55 + xi_57;
        const float xi_62 = -xi_26 * xi_61;
        const float xi_69 = u_0 * xi_68;
        const float xi_74 = u_0 * xi_73;
        const float xi_78 = xi_44 - xi_55 + xi_57;
        const float xi_85 = -xi_39;
        const float xi_96 = u_0 * xi_284;
        const float xi_97 = xi_96 * 0.25000000000000000f;
        const float xi_100 = xi_72 * xi_96;
        const float xi_119 = rho * (u_0 * u_0);
        const float xi_173 = rho * u_0;
        const float xi_174 = -vel0Term + xi_162 + xi_173 + xi_273 + xi_4;
        const float xi_175 = xi_140 * xi_174;
        const float xi_211 = xi_174 * xi_206;
        const float u_1 =
            xi_283 * xi_8 + xi_7 * (vel1Term + xi_17 + xi_20 + xi_276 + xi_9);
        const float xi_27 = u_1 * xi_283;
        const float xi_33 = xi_27 * 0.16666666666666667f;
        const float xi_46 = xi_27 * 0.083333333333333333f;
        const float xi_47 = omega_shear * xi_46;
        const float xi_48 = -xi_33 + xi_47;
        const float xi_63 = -xi_27 * xi_61;
        const float xi_70 = u_1 * 0.25000000000000000f;
        const float xi_71 = xi_271 * xi_70;
        const float xi_75 = u_1 * xi_72;
        const float xi_76 = xi_271 * xi_75;
        const float xi_77 = -xi_69 - xi_71 + xi_74 + xi_76;
        const float xi_79 = xi_69 + xi_71 - xi_74 - xi_76;
        const float xi_87 = xi_284 * xi_70;
        const float xi_89 = xi_284 * xi_75;
        const float xi_94 = -xi_46;
        const float xi_124 = rho * (u_1 * u_1);
        const float xi_125 = xi_10 + xi_123 + xi_124;
        const float xi_137 = rho * u_1;
        const float xi_139 =
            -vel1Term + xi_137 + xi_138 + xi_269 + xi_282 + xi_5;
        const float xi_141 = xi_139 * xi_140;
        const float xi_207 = xi_139 * xi_206;
        const float xi_217 = xi_216 * (u_0 * xi_137 + xi_138 + xi_273 + xi_9);
        const float xi_218 = -xi_215 - xi_217;
        const float xi_219 = xi_215 + xi_217;
        const float u_2 =
            xi_284 * xi_8 + xi_7 * (vel2Term + xi_22 + xi_25 + xi_274);
        const float xi_28 = u_2 * xi_284;
        const float xi_34 = xi_28 * 0.16666666666666667f;
        const float xi_35 = xi_28 * 0.083333333333333333f;
        const float xi_36 = omega_shear * xi_35;
        const float xi_37 = -xi_34 + xi_36;
        const float xi_42 =
            -omega_shear * xi_33 + xi_27 * 0.33333333333333333f + xi_37 + xi_41;
        const float xi_49 =
            -omega_shear * xi_38 + xi_26 * 0.33333333333333333f + xi_37 + xi_48;
        const float xi_53 =
            -omega_shear * xi_34 + xi_28 * 0.33333333333333333f + xi_41 + xi_48;
        const float xi_59 = -xi_35;
        const float xi_64 = -xi_28 * xi_61;
        const float xi_65 = -xi_27 * xi_56 + xi_33 + xi_62 + xi_63 + xi_64;
        const float xi_66 = xi_31 - xi_60 + xi_65;
        const float xi_67 = xi_36 + xi_59 + xi_66;
        const float xi_80 = -xi_31 + xi_60 + xi_65;
        const float xi_81 = xi_36 + xi_59 + xi_80;
        const float xi_83 = -xi_28 * xi_56 + xi_34;
        const float xi_84 = xi_51 - xi_82 + xi_83;
        const float xi_86 = xi_40 + xi_66 + xi_85;
        const float xi_88 = u_2 * xi_68;
        const float xi_90 = u_2 * xi_73;
        const float xi_91 = xi_87 + xi_88 - xi_89 - xi_90;
        const float xi_92 = xi_40 + xi_80 + xi_85;
        const float xi_93 = -xi_87 - xi_88 + xi_89 + xi_90;
        const float xi_95 = xi_47 + xi_62 + xi_63 + xi_64 + xi_84 + xi_94;
        const float xi_98 = u_2 * xi_271;
        const float xi_99 = xi_98 * 0.25000000000000000f;
        const float xi_101 = xi_72 * xi_98;
        const float xi_102 = xi_100 + xi_101 - xi_97 - xi_99;
        const float xi_103 = -xi_100 - xi_101 + xi_97 + xi_99;
        const float xi_104 = -xi_51 + xi_82 + xi_83;
        const float xi_105 = xi_104 + xi_47 + xi_62 + xi_63 + xi_64 + xi_94;
        const float xi_120 = rho * (u_2 * u_2);
        const float xi_128 = omega_bulk * (xi_119 + xi_120 + xi_122 + xi_125 +
                                           xi_127 + xi_18 + xi_23 + xi_278);
        const float xi_163 = -xi_120 + xi_286 + xi_289;
        const float xi_164 =
            omega_shear * (xi_0 + xi_125 + xi_162 + xi_163 + xi_17 - xi_268);
        const float xi_165 = xi_164 * 0.12500000000000000f;
        const float xi_167 =
            omega_shear *
            (xi_10 + xi_116 + xi_119 * 2.0f + xi_123 - xi_124 + xi_127 +
             xi_163 + xi_268 + xi_270 * -2.0f + xi_280 * -2.0f + xi_285);
        const float xi_169 =
            xi_167 * -0.041666666666666667f + xi_168 * -0.16666666666666667f;
        const float xi_170 = xi_109 * -0.10000000000000000f +
                             xi_115 * -0.050000000000000000f + xi_169;
        const float xi_171 = xi_112 * 0.028571428571428571f +
                             xi_118 * 0.014285714285714286f + xi_161 + xi_165 +
                             xi_166 + xi_170;
        const float xi_186 = xi_112 * -0.071428571428571429f +
                             xi_118 * -0.035714285714285714f + xi_166 +
                             xi_167 * 0.083333333333333333f +
                             xi_168 * 0.33333333333333333f;
        const float xi_191 =
            rho * u_2 - vel2Term + xi_121 + xi_126 + xi_188 + xi_279 + xi_6;
        const float xi_192 = xi_140 * xi_191;
        const float xi_200 =
            xi_108 * 0.095238095238095238f + xi_112 * -0.042857142857142857f +
            xi_118 * -0.021428571428571429f + xi_132 * 0.015873015873015873f -
            xi_161 - xi_165 + xi_170;
        const float xi_203 = xi_164 * 0.062500000000000000f;
        const float xi_208 =
            xi_111 * 0.083333333333333333f + xi_128 * 0.041666666666666667f;
        const float xi_209 = xi_207 + xi_208;
        const float xi_210 =
            xi_172 + xi_202 + xi_203 + xi_204 + xi_205 + xi_209;
        const float xi_212 =
            xi_167 * 0.020833333333333333f + xi_168 * 0.083333333333333333f;
        const float xi_213 = -xi_211 + xi_212;
        const float xi_214 = xi_187 + xi_213;
        const float xi_220 = xi_211 + xi_212;
        const float xi_221 = xi_185 + xi_220;
        const float xi_222 = -xi_207 + xi_208;
        const float xi_223 =
            xi_157 + xi_202 + xi_203 + xi_204 + xi_205 + xi_222;
        const float xi_226 = xi_216 * (u_2 * xi_137 + xi_133 + xi_18 + xi_281);
        const float xi_230 =
            xi_169 + xi_224 + xi_225 + xi_226 + xi_227 + xi_228 + xi_229;
        const float xi_239 = xi_191 * xi_206;
        const float xi_241 = xi_239 + xi_240;
        const float xi_242 = -xi_232 + xi_234 - xi_236 + xi_238 + xi_241;
        const float xi_247 = xi_209 - xi_243 + xi_244 - xi_245 + xi_246;
        const float xi_248 = xi_222 + xi_243 - xi_244 + xi_245 - xi_246;
        const float xi_249 =
            xi_169 - xi_224 + xi_225 - xi_226 + xi_227 + xi_228 + xi_229;
        const float xi_251 = -xi_203;
        const float xi_254 =
            xi_201 + xi_208 + xi_241 + xi_250 + xi_251 + xi_252 + xi_253;
        const float xi_256 = xi_216 * (u_2 * xi_173 + xi_11 + xi_176 + xi_279);
        const float xi_257 = -xi_255 - xi_256;
        const float xi_262 = xi_213 - xi_258 + xi_259 - xi_260 + xi_261;
        const float xi_263 = xi_255 + xi_256;
        const float xi_264 = xi_220 + xi_258 - xi_259 + xi_260 - xi_261;
        const float xi_265 = -xi_239 + xi_240;
        const float xi_266 = xi_232 - xi_234 + xi_236 - xi_238 + xi_265;
        const float xi_267 =
            xi_199 + xi_208 + xi_250 + xi_251 + xi_252 + xi_253 + xi_265;
        const float forceTerm_0 = xi_26 * xi_29 - xi_26 + xi_27 * xi_29 -
                                  xi_27 + xi_28 * xi_29 - xi_28;
        const float forceTerm_1 = xi_30 - xi_32 + xi_42;
        const float forceTerm_2 = -xi_30 + xi_32 + xi_42;
        const float forceTerm_3 = -xi_43 + xi_45 + xi_49;
        const float forceTerm_4 = xi_43 - xi_45 + xi_49;
        const float forceTerm_5 = xi_50 - xi_52 + xi_53;
        const float forceTerm_6 = -xi_50 + xi_52 + xi_53;
        const float forceTerm_7 = xi_58 + xi_67 + xi_77;
        const float forceTerm_8 = xi_67 + xi_78 + xi_79;
        const float forceTerm_9 = xi_58 + xi_79 + xi_81;
        const float forceTerm_10 = xi_77 + xi_78 + xi_81;
        const float forceTerm_11 = xi_84 + xi_86 + xi_91;
        const float forceTerm_12 = xi_84 + xi_92 + xi_93;
        const float forceTerm_13 = xi_102 + xi_58 + xi_95;
        const float forceTerm_14 = xi_103 + xi_78 + xi_95;
        const float forceTerm_15 = xi_104 + xi_86 + xi_93;
        const float forceTerm_16 = xi_104 + xi_91 + xi_92;
        const float forceTerm_17 = xi_103 + xi_105 + xi_58;
        const float forceTerm_18 = xi_102 + xi_105 + xi_78;
        _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_0 + xi_108 * 0.14285714285714286f +
            xi_109 * 0.20000000000000000f - xi_111 +
            xi_112 * 0.085714285714285714f + xi_115 * 0.10000000000000000f +
            xi_118 * 0.042857142857142857f + xi_128 * -0.50000000000000000f +
            xi_132 * 0.023809523809523810f + xi_278;
        _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_1 - xi_136 + xi_141 - xi_146 + xi_157 + xi_171 + xi_268;
        _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_2 + xi_136 - xi_141 + xi_146 + xi_171 + xi_172 + xi_285;
        _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_3 - xi_175 + xi_178 + xi_180 + xi_185 + xi_186 + xi_280;
        _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_4 + xi_175 - xi_178 - xi_180 + xi_186 + xi_187 + xi_270;
        _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_5 - xi_190 + xi_192 - xi_194 + xi_199 + xi_200 + xi_286;
        _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_6 + xi_190 - xi_192 + xi_194 + xi_200 + xi_201 + xi_289;
        _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_7 + xi_210 + xi_214 + xi_218 + xi_273;
        _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_8 + xi_210 + xi_219 + xi_221 + xi_276;
        _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_9 + xi_214 + xi_219 + xi_223 + xi_269;
        _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_10 + xi_218 + xi_221 + xi_223 + xi_288;
        _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_11 + xi_230 + xi_242 + xi_247 + xi_275;
        _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_12 + xi_242 + xi_248 + xi_249 + xi_282;
        _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_13 + xi_254 + xi_257 + xi_262 + xi_287;
        _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_14 + xi_254 + xi_263 + xi_264 + xi_274;
        _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_15 + xi_247 + xi_249 + xi_266 + xi_281;
        _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_16 + xi_230 + xi_248 + xi_266 + xi_272;
        _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_17 + xi_262 + xi_263 + xi_267 + xi_277;
        _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_18 + xi_257 + xi_264 + xi_267 + xi_279;
      }
    }
  }
}
} // namespace internal_69764eed2d0964e29e3b97d1054b4693

void CollideSweepSinglePrecisionThermalized::run(IBlock *block) {
  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);
  auto force = block->getData<field::GhostLayerField<float, 3>>(forceID);

  auto block_offset_0 = this->block_offset_0_;
  auto &kT = this->kT_;
  auto &omega_shear = this->omega_shear_;
  auto block_offset_2 = this->block_offset_2_;
  auto &omega_bulk = this->omega_bulk_;
  auto &omega_odd = this->omega_odd_;
  auto &omega_even = this->omega_even_;
  auto &time_step = this->time_step_;
  auto block_offset_1 = this->block_offset_1_;
  auto &seed = this->seed_;
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
  internal_69764eed2d0964e29e3b97d1054b4693::
      collidesweepsingleprecisionthermalized_collidesweepsingleprecisionthermalized(
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

  auto block_offset_0 = this->block_offset_0_;
  auto &kT = this->kT_;
  auto &omega_shear = this->omega_shear_;
  auto block_offset_2 = this->block_offset_2_;
  auto &omega_bulk = this->omega_bulk_;
  auto &omega_odd = this->omega_odd_;
  auto &omega_even = this->omega_even_;
  auto &time_step = this->time_step_;
  auto block_offset_1 = this->block_offset_1_;
  auto &seed = this->seed_;
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
  internal_69764eed2d0964e29e3b97d1054b4693::
      collidesweepsingleprecisionthermalized_collidesweepsingleprecisionthermalized(
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