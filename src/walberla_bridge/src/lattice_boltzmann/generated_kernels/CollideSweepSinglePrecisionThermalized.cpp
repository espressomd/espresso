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

// kernel generated with pystencils v1.3.1+2.g60e24c4, lbmpy v1.3.1+6.gcd1bc2f.dirty, lbmpy_walberla/pystencils_walberla from waLBerla commit 065ce5f311850371a97ac4766f47dbb5ca8424ba

#include <cmath>

#include "CollideSweepSinglePrecisionThermalized.h"
#include "core/DataTypes.h"
#include "core/Macros.h"

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

namespace internal_69764eed2d0964e29e3b97d1054b4693 {
static FUNC_PREFIX void collidesweepsingleprecisionthermalized_collidesweepsingleprecisionthermalized(float *RESTRICT const _data_force, float *RESTRICT _data_pdfs, int64_t const _size_force_0, int64_t const _size_force_1, int64_t const _size_force_2, int64_t const _stride_force_0, int64_t const _stride_force_1, int64_t const _stride_force_2, int64_t const _stride_force_3, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, uint32_t block_offset_0, uint32_t block_offset_1, uint32_t block_offset_2, float kT, float omega_bulk, float omega_even, float omega_odd, float omega_shear, uint32_t seed, uint32_t time_step) {
  const float xi_21 = omega_bulk * 0.5f;
  const float xi_48 = omega_shear * 0.041666666666666664f;
  const float xi_52 = omega_bulk * 0.041666666666666664f;
  const float xi_63 = omega_shear * 0.125f;
  const float xi_98 = 3.7416573867739413f;
  const float xi_101 = 5.4772255750516612f;
  const float xi_105 = 2.4494897427831779f;
  const float xi_108 = 8.3666002653407556f;
  const float xi_151 = omega_odd * 0.25f;
  const float xi_161 = omega_odd * 0.083333333333333329f;
  const float xi_174 = 1.7320508075688772f;
  const float xi_218 = omega_shear * 0.25f;
  const float xi_224 = omega_odd * 0.041666666666666664f;
  const float xi_227 = omega_odd * 0.125f;
  const float rr_0 = 0.0f;
  const float xi_46 = rr_0 * 0.041666666666666664f;
  for (int64_t ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1) {
    float *RESTRICT _data_pdfs_20_37 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_311 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_39 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_34 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_31 = _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_310 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_38 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_317 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_32 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_32 = _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    float *RESTRICT _data_pdfs_20_316 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    float *RESTRICT _data_pdfs_20_36 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_312 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_315 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_31 = _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    float *RESTRICT _data_pdfs_20_35 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_313 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_33 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_314 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    float *RESTRICT _data_pdfs_20_318 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1) {
      float *RESTRICT _data_pdfs_20_37_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_37;
      float *RESTRICT _data_pdfs_20_311_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_311;
      float *RESTRICT _data_pdfs_20_39_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_39;
      float *RESTRICT _data_pdfs_20_34_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      float *RESTRICT _data_pdfs_20_31_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_31;
      float *RESTRICT _data_pdfs_20_310_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_310;
      float *RESTRICT _data_pdfs_20_38_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_38;
      float *RESTRICT _data_pdfs_20_317_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_317;
      float *RESTRICT _data_pdfs_20_32_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_32;
      float *RESTRICT _data_force_20_32_10 = _stride_force_1 * ctr_1 + _data_force_20_32;
      float *RESTRICT _data_pdfs_20_316_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_316;
      float *RESTRICT _data_pdfs_20_30_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      float *RESTRICT _data_pdfs_20_36_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_36;
      float *RESTRICT _data_pdfs_20_312_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_312;
      float *RESTRICT _data_pdfs_20_315_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_315;
      float *RESTRICT _data_force_20_31_10 = _stride_force_1 * ctr_1 + _data_force_20_31;
      float *RESTRICT _data_pdfs_20_35_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_35;
      float *RESTRICT _data_pdfs_20_313_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_313;
      float *RESTRICT _data_pdfs_20_33_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      float *RESTRICT _data_pdfs_20_314_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_314;
      float *RESTRICT _data_force_20_30_10 = _stride_force_1 * ctr_1 + _data_force_20_30;
      float *RESTRICT _data_pdfs_20_318_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_318;
      for (int64_t ctr_0 = 0; ctr_0 < _size_force_0; ctr_0 += 1) {
        const float xi_264 = _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0];
        const float xi_265 = _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0];
        const float xi_266 = _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0];
        const float xi_267 = _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0];
        const float xi_268 = _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0];
        const float xi_269 = _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0];
        const float xi_270 = _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0];
        const float xi_271 = _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0];
        const float xi_272 = _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0];
        const float xi_273 = _data_force_20_32_10[_stride_force_0 * ctr_0];
        const float xi_274 = _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0];
        const float xi_275 = _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0];
        const float xi_276 = _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0];
        const float xi_277 = _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0];
        const float xi_278 = _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0];
        const float xi_279 = _data_force_20_31_10[_stride_force_0 * ctr_0];
        const float xi_280 = _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0];
        const float xi_281 = _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0];
        const float xi_282 = _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0];
        const float xi_283 = _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0];
        const float xi_284 = _data_force_20_30_10[_stride_force_0 * ctr_0];
        const float xi_285 = _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0];

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
        const float xi_2 = xi_269 + xi_270;
        const float xi_3 = xi_2 + xi_283 + xi_285;
        const float xi_4 = xi_264 + xi_268;
        const float xi_5 = xi_265 + xi_278;
        const float xi_6 = xi_277 + xi_280;
        const float xi_7 = xi_271 + xi_282;
        const float xi_8 = xi_266 + xi_7;
        const float xi_9 = xi_272 + xi_274;
        const float xi_12 = xi_264 + xi_281;
        const float xi_13 = xi_269 + xi_270 * -1.0f;
        const float xi_14 = xi_266 + xi_277;
        const float xi_15 = xi_265 * -1.0f;
        const float xi_16 = xi_283 * -1.0f;
        const float xi_17 = xi_271 + xi_274;
        const float xi_22 = xi_279 * 0.16666666666666666f;
        const float xi_23 = xi_279 * 0.083333333333333329f;
        const float xi_34 = xi_284 * 0.16666666666666666f;
        const float xi_35 = xi_284 * 0.083333333333333329f;
        const float xi_40 = xi_273 * 0.16666666666666666f;
        const float xi_41 = xi_273 * 0.083333333333333329f;
        const float xi_59 = xi_279 * 0.25f;
        const float xi_64 = xi_279 * xi_63;
        const float xi_99 = random_3_0 - 0.5f;
        const float xi_102 = random_3_2 - 0.5f;
        const float xi_104 = random_1_1 - 0.5f;
        const float xi_109 = random_3_1 - 0.5f;
        const float xi_113 = xi_275 * -1.0f;
        const float xi_119 = xi_12 + xi_3;
        const float xi_122 = xi_265 * 2.0f + xi_267 * 5.0f + xi_274 * 2.0f + xi_277 * 2.0f + xi_278 * 2.0f + xi_282 * 5.0f;
        const float xi_125 = xi_281 * 2.0f;
        const float xi_126 = xi_283 * 2.0f;
        const float xi_127 = xi_271 * 2.0f + xi_285 * 2.0f;
        const float xi_130 = xi_268 * -1.0f;
        const float xi_134 = random_0_1 - 0.5f;
        const float xi_138 = xi_271 * -1.0f;
        const float xi_139 = xi_285 * -1.0f;
        const float xi_140 = xi_280 * -1.0f;
        const float xi_141 = xi_276 * -1.0f;
        const float xi_146 = xi_278 * -1.0f;
        const float xi_147 = xi_15 + xi_277;
        const float xi_148 = xi_272 * -1.0f;
        const float xi_149 = xi_148 + xi_268;
        const float xi_150 = xi_146 + xi_147 + xi_149 + xi_274;
        const float xi_152 = random_2_3 - 0.5f;
        const float xi_156 = xi_274 * -1.0f;
        const float xi_157 = xi_277 * -1.0f;
        const float xi_158 = xi_270 * 2.0f;
        const float xi_159 = xi_264 * -2.0f + xi_269 * 2.0f;
        const float xi_160 = xi_149 + xi_156 + xi_157 + xi_158 * -1.0f + xi_159 + xi_266 * 2.0f + xi_5;
        const float xi_162 = xi_160 * xi_161;
        const float xi_163 = random_1_2 - 0.5f;
        const float xi_175 = random_0_0 - 0.5f;
        const float xi_186 = xi_16 + xi_281;
        const float xi_187 = xi_139 + xi_186;
        const float xi_188 = xi_187 + xi_267 + xi_271 + xi_282 * -1.0f;
        const float xi_189 = random_2_1 - 0.5f;
        const float xi_190 = xi_267 * -1.0f;
        const float xi_191 = xi_158 * -1.0f + xi_159 * -1.0f + xi_187 * -1.0f + xi_190 * -1.0f + xi_266 * 2.0f + xi_7 * -1.0f;
        const float xi_192 = xi_161 * xi_191;
        const float xi_193 = random_2_0 - 0.5f;
        const float xi_199 = xi_146 + xi_156 + xi_265;
        const float xi_200 = xi_140 * -1.0f + xi_199 * -1.0f + xi_276 * -1.0f + xi_277 * -1.0f;
        const float xi_201 = random_2_2 - 0.5f;
        const float xi_202 = xi_125 * -1.0f + xi_126 * -1.0f + xi_127 + xi_141 + xi_199 + xi_6;
        const float xi_203 = xi_161 * xi_202;
        const float xi_204 = random_1_3 - 0.5f;
        const float xi_225 = xi_160 * xi_224;
        const float xi_228 = xi_150 * xi_227;
        const float xi_234 = xi_202 * xi_224;
        const float xi_235 = xi_200 * xi_227;
        const float xi_255 = xi_188 * xi_227;
        const float xi_256 = xi_191 * xi_224;
        const float xi_24 = rr_0 * xi_23;
        const float xi_36 = rr_0 * xi_35;
        const float xi_42 = rr_0 * xi_41;
        const float xi_47 = xi_284 * xi_46;
        const float xi_51 = xi_279 * xi_46;
        const float xi_73 = xi_273 * xi_46;
        const float vel0Term = xi_267 + xi_3;
        const float vel1Term = xi_4 + xi_5;
        const float vel2Term = xi_281 + xi_6;
        const float delta_rho = vel0Term + vel1Term + vel2Term + xi_275 + xi_276 + xi_8 + xi_9;
        const float xi_10 = delta_rho + 1.0f;
        const float xi_96 = kT * xi_10;
        const float xi_97 = powf(xi_96 * (-1.0f * ((omega_even * -1.0f + 1.0f) * (omega_even * -1.0f + 1.0f)) + 1.0f), 0.5f);
        const float xi_100 = xi_97 * xi_98 * xi_99;
        const float xi_103 = xi_101 * xi_102 * xi_97;
        const float xi_106 = powf(xi_96 * (-1.0f * ((omega_bulk * -1.0f + 1.0f) * (omega_bulk * -1.0f + 1.0f)) + 1.0f), 0.5f);
        const float xi_107 = xi_104 * xi_105 * xi_106;
        const float xi_110 = xi_108 * xi_109 * xi_97;
        const float xi_132 = xi_100 * 0.11904761904761904f;
        const float xi_135 = powf(xi_96 * (-1.0f * ((omega_shear * -1.0f + 1.0f) * (omega_shear * -1.0f + 1.0f)) + 1.0f), 0.5f);
        const float xi_136 = xi_135 * 0.5f;
        const float xi_137 = xi_134 * xi_136;
        const float xi_153 = powf(xi_96 * (-1.0f * ((omega_odd * -1.0f + 1.0f) * (omega_odd * -1.0f + 1.0f)) + 1.0f), 0.5f);
        const float xi_154 = xi_153 * 1.4142135623730951f;
        const float xi_155 = xi_154 * 0.5f;
        const float xi_164 = xi_105 * xi_153;
        const float xi_165 = xi_164 * 0.16666666666666666f;
        const float xi_166 = xi_163 * xi_165;
        const float xi_167 = xi_162 + xi_166;
        const float xi_168 = xi_150 * xi_151 + xi_152 * xi_155 + xi_167;
        const float xi_170 = xi_103 * 0.10000000000000001f;
        const float xi_176 = xi_135 * xi_174 * xi_175;
        const float xi_177 = xi_176 * 0.16666666666666666f;
        const float xi_185 = xi_110 * 0.071428571428571425f;
        const float xi_194 = xi_165 * xi_193;
        const float xi_195 = xi_192 + xi_194;
        const float xi_196 = xi_151 * xi_188 + xi_155 * xi_189 + xi_195;
        const float xi_198 = xi_110 * 0.042857142857142858f;
        const float xi_205 = xi_165 * xi_204;
        const float xi_206 = xi_203 + xi_205;
        const float xi_207 = xi_151 * xi_200 + xi_155 * xi_201 + xi_206;
        const float xi_208 = xi_134 * xi_135 * 0.25f;
        const float xi_211 = xi_100 * 0.083333333333333329f;
        const float xi_215 = xi_192 * -1.0f + xi_194 * -1.0f;
        const float xi_216 = xi_136 * (random_0_2 - 0.5f);
        const float xi_223 = xi_136 * (random_1_0 - 0.5f);
        const float xi_229 = xi_164 * 0.083333333333333329f;
        const float xi_230 = xi_163 * xi_229;
        const float xi_231 = xi_154 * 0.25f;
        const float xi_232 = xi_152 * xi_231;
        const float xi_236 = xi_204 * xi_229;
        const float xi_237 = xi_201 * xi_231;
        const float xi_238 = xi_234 * -1.0f + xi_235 + xi_236 * -1.0f + xi_237;
        const float xi_240 = xi_110 * 0.014285714285714285f;
        const float xi_242 = xi_100 * 0.023809523809523808f;
        const float xi_245 = xi_234 + xi_235 * -1.0f + xi_236 + xi_237 * -1.0f;
        const float xi_247 = xi_208 * -1.0f;
        const float xi_250 = xi_110 * 0.035714285714285712f;
        const float xi_252 = xi_136 * (random_0_3 - 0.5f);
        const float xi_257 = xi_189 * xi_231;
        const float xi_258 = xi_193 * xi_229;
        const float xi_259 = xi_255 * -1.0f + xi_256 + xi_257 * -1.0f + xi_258;
        const float xi_261 = xi_255 + xi_256 * -1.0f + xi_257 + xi_258 * -1.0f;
        const float rho = xi_10;
        const float xi_0 = ((1.0f) / (rho));
        const float xi_11 = xi_0 * 0.5f;
        const float u_0 = xi_0 * (vel0Term + xi_12 * -1.0f + xi_8 * -1.0f) + xi_11 * xi_284;
        const float xi_18 = u_0 * xi_284;
        const float xi_29 = xi_18 * 0.16666666666666666f;
        const float xi_30 = xi_29 * -1.0f;
        const float xi_31 = xi_18 * 0.083333333333333329f;
        const float xi_32 = omega_shear * xi_31 + xi_30;
        const float xi_49 = xi_18 * xi_48 + xi_30;
        const float xi_50 = xi_35 + xi_47 * -1.0f + xi_49;
        const float xi_53 = xi_18 * xi_52;
        const float xi_60 = u_0 * xi_59;
        const float xi_65 = u_0 * xi_64;
        const float xi_69 = xi_35 * -1.0f + xi_47 + xi_49;
        const float xi_76 = omega_shear * u_0 * xi_284 * -0.083333333333333329f;
        const float xi_86 = u_0 * xi_273;
        const float xi_87 = xi_86 * 0.25f;
        const float xi_90 = xi_63 * xi_86;
        const float xi_112 = u_0 * u_0;
        const float u_1 = xi_0 * (vel1Term + xi_13 * -1.0f + xi_14 * -1.0f + xi_9 * -1.0f) + xi_11 * xi_279;
        const float xi_19 = u_1 * xi_279;
        const float xi_27 = xi_19 * 0.16666666666666666f;
        const float xi_37 = omega_shear * u_1 * xi_279 * -0.083333333333333329f;
        const float xi_43 = xi_27 * -1.0f;
        const float xi_44 = xi_19 * 0.083333333333333329f;
        const float xi_54 = xi_19 * xi_52;
        const float xi_61 = u_1 * 0.25f;
        const float xi_62 = xi_284 * xi_61;
        const float xi_66 = u_1 * xi_63;
        const float xi_67 = xi_284 * xi_66;
        const float xi_68 = xi_60 + xi_62 + xi_65 * -1.0f + xi_67 * -1.0f;
        const float xi_70 = xi_60 * -1.0f + xi_62 * -1.0f + xi_65 + xi_67;
        const float xi_78 = xi_273 * xi_61;
        const float xi_80 = xi_273 * xi_66;
        const float xi_111 = rho * (u_1 * u_1);
        const float xi_118 = xi_111 * -1.0f;
        const float xi_217 = rho * u_1;
        const float xi_219 = xi_218 * (u_0 * xi_217 + xi_13 + xi_264 + xi_266 * -1.0f);
        const float xi_220 = xi_216 * -1.0f + xi_219 * -1.0f;
        const float xi_221 = xi_216 + xi_219;
        const float u_2 = xi_0 * (vel2Term + xi_15 * -1.0f + xi_16 * -1.0f + xi_17 * -1.0f + xi_276 * -1.0f + xi_278 * -1.0f + xi_285 * -1.0f) + xi_11 * xi_273;
        const float xi_20 = u_2 * xi_273;
        const float xi_25 = xi_20 * 0.16666666666666666f;
        const float xi_26 = xi_25 * -1.0f;
        const float xi_28 = xi_20 * 0.083333333333333329f;
        const float xi_33 = omega_shear * xi_27 * -1.0f + omega_shear * xi_28 + xi_19 * 0.33333333333333331f + xi_26 + xi_32;
        const float xi_38 = omega_shear * u_2 * xi_273 * -0.083333333333333329f;
        const float xi_39 = omega_shear * xi_29 + u_0 * xi_284 * -0.33333333333333331f + xi_25 + xi_27 + xi_37 + xi_38;
        const float xi_45 = omega_shear * xi_25 * -1.0f + omega_shear * xi_44 + xi_20 * 0.33333333333333331f + xi_32 + xi_43;
        const float xi_55 = xi_20 * xi_52;
        const float xi_56 = xi_19 * xi_48 + xi_43 + xi_53 + xi_54 + xi_55;
        const float xi_57 = xi_23 * -1.0f + xi_51 + xi_56;
        const float xi_58 = xi_28 + xi_38 + xi_57;
        const float xi_71 = xi_23 + xi_51 * -1.0f + xi_56;
        const float xi_72 = xi_28 + xi_38 + xi_71;
        const float xi_74 = xi_20 * xi_48 + xi_26;
        const float xi_75 = xi_41 * -1.0f + xi_73 + xi_74;
        const float xi_77 = xi_31 + xi_57 + xi_76;
        const float xi_79 = u_2 * xi_59;
        const float xi_81 = u_2 * xi_64;
        const float xi_82 = xi_78 * -1.0f + xi_79 * -1.0f + xi_80 + xi_81;
        const float xi_83 = xi_31 + xi_71 + xi_76;
        const float xi_84 = xi_78 + xi_79 + xi_80 * -1.0f + xi_81 * -1.0f;
        const float xi_85 = xi_37 + xi_44 + xi_53 + xi_54 + xi_55 + xi_75;
        const float xi_88 = u_2 * xi_284;
        const float xi_89 = xi_88 * 0.25f;
        const float xi_91 = xi_63 * xi_88;
        const float xi_92 = xi_87 + xi_89 + xi_90 * -1.0f + xi_91 * -1.0f;
        const float xi_93 = xi_87 * -1.0f + xi_89 * -1.0f + xi_90 + xi_91;
        const float xi_94 = xi_41 + xi_73 * -1.0f + xi_74;
        const float xi_95 = xi_37 + xi_44 + xi_53 + xi_54 + xi_55 + xi_94;
        const float xi_114 = rho * (u_2 * u_2);
        const float xi_115 = xi_113 + xi_114 * 0.66666666666666663f + xi_276 * 3.0f + xi_280 * 3.0f;
        const float xi_116 = rho * xi_112 * 1.6666666666666667f + xi_111 * 0.66666666666666663f + xi_115 + xi_265 * -3.0f + xi_268 * 3.0f + xi_272 * 3.0f + xi_274 * -3.0f + xi_277 * -3.0f + xi_278 * -3.0f;
        const float xi_117 = omega_even * xi_116;
        const float xi_120 = rho * xi_112 + xi_113 * -1.0f + xi_114 + xi_118 * -1.0f + xi_119 * -1.0f + xi_14 * -1.0f + xi_17 * -1.0f + xi_5 * -1.0f;
        const float xi_121 = omega_bulk * xi_120;
        const float xi_123 = xi_111 * 2.3333333333333335f + xi_115 + xi_122 + xi_268 * -2.0f + xi_271 * -5.0f + xi_272 * -2.0f + xi_281 * -5.0f + xi_283 * -5.0f + xi_285 * -5.0f;
        const float xi_124 = omega_even * xi_123;
        const float xi_128 = xi_113 + xi_114 * 3.0f + xi_122 + xi_125 + xi_126 + xi_127 + xi_264 * -7.0f + xi_266 * -7.0f + xi_268 * 5.0f + xi_269 * -7.0f + xi_270 * -7.0f + xi_272 * 5.0f + xi_276 * -4.0f + xi_280 * -4.0f;
        const float xi_129 = omega_even * xi_128;
        const float xi_131 = xi_129 * 0.01984126984126984f;
        const float xi_133 = xi_131 + xi_132;
        const float xi_142 = xi_114 + xi_140 + xi_141 + xi_266;
        const float xi_143 = omega_shear * (xi_118 * -1.0f + xi_138 * -1.0f + xi_139 * -1.0f + xi_142 * -1.0f + xi_16 * -1.0f + xi_2 * -1.0f + xi_272 * -1.0f + xi_281 + xi_4 * -1.0f);
        const float xi_144 = xi_143 * 0.125f;
        const float xi_145 = xi_137 * -1.0f + xi_144 * -1.0f;
        const float xi_169 = xi_117 * 0.050000000000000003f;
        const float xi_171 = rho * xi_112 * 2.0f + xi_111 * -1.0f + xi_119 * -1.0f + xi_130 * -1.0f + xi_142 * -1.0f + xi_148 * -1.0f + xi_265 * 2.0f + xi_267 * -2.0f + xi_271 * -1.0f + xi_274 * 2.0f + xi_277 * 2.0f + xi_278 * 2.0f + xi_282 * -2.0f;
        const float xi_172 = omega_shear * xi_171;
        const float xi_173 = xi_172 * 0.041666666666666664f;
        const float xi_178 = xi_173 + xi_177;
        const float xi_179 = xi_169 + xi_170 + xi_178;
        const float xi_180 = xi_131 * -1.0f + xi_132 * -1.0f;
        const float xi_181 = xi_137 + xi_144;
        const float xi_182 = xi_173 * -1.0f + xi_177 * -1.0f;
        const float xi_183 = xi_169 * -1.0f + xi_170 * -1.0f + xi_182;
        const float xi_184 = xi_124 * 0.035714285714285712f;
        const float xi_197 = xi_124 * 0.021428571428571429f;
        const float xi_209 = xi_143 * 0.0625f;
        const float xi_210 = xi_129 * 0.013888888888888888f;
        const float xi_212 = xi_107 * 0.083333333333333329f + xi_121 * 0.041666666666666664f;
        const float xi_213 = xi_172 * 0.020833333333333332f + xi_176 * 0.083333333333333329f + xi_212;
        const float xi_214 = xi_167 + xi_208 + xi_209 + xi_210 + xi_211 + xi_213;
        const float xi_222 = xi_162 * -1.0f + xi_166 * -1.0f + xi_208 + xi_209 + xi_210 + xi_211 + xi_213;
        const float xi_226 = xi_218 * (u_2 * xi_217 + xi_147 + xi_156 + xi_278);
        const float xi_233 = xi_223 + xi_225 * -1.0f + xi_226 + xi_228 + xi_230 * -1.0f + xi_232;
        const float xi_239 = xi_124 * 0.0071428571428571426f;
        const float xi_241 = xi_129 * 0.003968253968253968f;
        const float xi_243 = xi_241 * -1.0f + xi_242 * -1.0f;
        const float xi_244 = xi_103 * 0.050000000000000003f + xi_117 * 0.025000000000000001f + xi_182 + xi_212 + xi_239 * -1.0f + xi_240 * -1.0f + xi_243;
        const float xi_246 = omega_bulk * xi_120 * -0.041666666666666664f + omega_even * xi_116 * -0.025000000000000001f + xi_101 * xi_102 * xi_97 * -0.050000000000000003f + xi_104 * xi_105 * xi_106 * -0.083333333333333329f + xi_178 + xi_239 + xi_240 + xi_241 + xi_242;
        const float xi_248 = xi_209 * -1.0f;
        const float xi_249 = xi_124 * 0.017857142857142856f;
        const float xi_251 = xi_206 + xi_213 + xi_243 + xi_247 + xi_248 + xi_249 + xi_250;
        const float xi_253 = xi_218 * (rho * u_0 * u_2 + xi_138 + xi_186 + xi_285);
        const float xi_254 = xi_252 * -1.0f + xi_253 * -1.0f;
        const float xi_260 = xi_252 + xi_253;
        const float xi_262 = xi_223 + xi_225 + xi_226 + xi_228 * -1.0f + xi_230 + xi_232 * -1.0f;
        const float xi_263 = xi_203 * -1.0f + xi_205 * -1.0f + xi_213 + xi_243 + xi_247 + xi_248 + xi_249 + xi_250;
        const float forceTerm_0 = xi_18 * xi_21 + xi_18 * -1.0f + xi_19 * xi_21 + xi_19 * -1.0f + xi_20 * xi_21 + xi_20 * -1.0f;
        const float forceTerm_1 = xi_22 + xi_24 * -1.0f + xi_33;
        const float forceTerm_2 = xi_22 * -1.0f + xi_24 + xi_33;
        const float forceTerm_3 = xi_34 * -1.0f + xi_36 + xi_39 * -1.0f;
        const float forceTerm_4 = xi_34 + xi_36 * -1.0f + xi_39 * -1.0f;
        const float forceTerm_5 = xi_40 + xi_42 * -1.0f + xi_45;
        const float forceTerm_6 = xi_40 * -1.0f + xi_42 + xi_45;
        const float forceTerm_7 = xi_50 * -1.0f + xi_58 * -1.0f + xi_68 * -1.0f;
        const float forceTerm_8 = xi_58 * -1.0f + xi_69 * -1.0f + xi_70 * -1.0f;
        const float forceTerm_9 = xi_50 * -1.0f + xi_70 * -1.0f + xi_72 * -1.0f;
        const float forceTerm_10 = xi_68 * -1.0f + xi_69 * -1.0f + xi_72 * -1.0f;
        const float forceTerm_11 = xi_75 * -1.0f + xi_77 * -1.0f + xi_82 * -1.0f;
        const float forceTerm_12 = xi_75 * -1.0f + xi_83 * -1.0f + xi_84 * -1.0f;
        const float forceTerm_13 = xi_50 * -1.0f + xi_85 * -1.0f + xi_92 * -1.0f;
        const float forceTerm_14 = xi_69 * -1.0f + xi_85 * -1.0f + xi_93 * -1.0f;
        const float forceTerm_15 = xi_77 * -1.0f + xi_84 * -1.0f + xi_94 * -1.0f;
        const float forceTerm_16 = xi_82 * -1.0f + xi_83 * -1.0f + xi_94 * -1.0f;
        const float forceTerm_17 = xi_50 * -1.0f + xi_93 * -1.0f + xi_95 * -1.0f;
        const float forceTerm_18 = xi_69 * -1.0f + xi_92 * -1.0f + xi_95 * -1.0f;
        _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0] = forceTerm_0 + xi_100 * 0.14285714285714285f + xi_103 * 0.20000000000000001f + xi_107 * -1.0f + xi_110 * 0.085714285714285715f + xi_117 * 0.10000000000000001f + xi_121 * -0.5f + xi_124 * 0.042857142857142858f + xi_129 * 0.023809523809523808f + xi_275;
        _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0] = forceTerm_1 + omega_even * xi_123 * 0.014285714285714285f + xi_108 * xi_109 * xi_97 * 0.028571428571428571f + xi_130 * -1.0f + xi_133 * -1.0f + xi_145 * -1.0f + xi_168 * -1.0f + xi_179 * -1.0f;
        _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0] = forceTerm_2 + xi_110 * 0.028571428571428571f + xi_124 * 0.014285714285714285f + xi_168 + xi_180 + xi_181 + xi_183 + xi_272;
        _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0] = forceTerm_3 + xi_172 * 0.083333333333333329f + xi_176 * 0.33333333333333331f + xi_180 + xi_184 * -1.0f + xi_185 * -1.0f + xi_196 + xi_282;
        _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0] = forceTerm_4 + omega_shear * xi_171 * 0.083333333333333329f + xi_133 * -1.0f + xi_135 * xi_174 * xi_175 * 0.33333333333333331f + xi_184 * -1.0f + xi_185 * -1.0f + xi_190 * -1.0f + xi_196 * -1.0f;
        _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0] = forceTerm_5 + omega_even * xi_128 * 0.015873015873015872f + xi_140 * -1.0f + xi_179 * -1.0f + xi_181 * -1.0f + xi_197 * -1.0f + xi_198 * -1.0f + xi_207 * -1.0f + xi_97 * xi_98 * xi_99 * 0.095238095238095233f;
        _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0] = forceTerm_6 + xi_100 * 0.095238095238095233f + xi_129 * 0.015873015873015872f + xi_145 + xi_183 + xi_197 * -1.0f + xi_198 * -1.0f + xi_207 + xi_276;
        _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0] = forceTerm_7 + xi_214 + xi_215 + xi_220 + xi_264;
        _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0] = forceTerm_8 + xi_195 + xi_214 + xi_221 + xi_270;
        _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0] = forceTerm_9 + xi_215 + xi_221 + xi_222 + xi_266;
        _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0] = forceTerm_10 + xi_195 + xi_220 + xi_222 + xi_269;
        _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0] = forceTerm_11 + xi_233 + xi_238 + xi_244 + xi_265;
        _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0] = forceTerm_12 + xi_157 * -1.0f + xi_233 * -1.0f + xi_245 * -1.0f + xi_246 * -1.0f;
        _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0] = forceTerm_13 + xi_251 + xi_254 + xi_259 + xi_281;
        _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0] = forceTerm_14 + xi_251 + xi_260 + xi_261 + xi_283;
        _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0] = forceTerm_15 + xi_146 * -1.0f + xi_238 * -1.0f + xi_246 * -1.0f + xi_262 * -1.0f;
        _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0] = forceTerm_16 + xi_244 + xi_245 + xi_262 + xi_274;
        _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0] = forceTerm_17 + xi_259 + xi_260 + xi_263 + xi_271;
        _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0] = forceTerm_18 + xi_254 + xi_261 + xi_263 + xi_285;
      }
    }
  }
}
} // namespace internal_69764eed2d0964e29e3b97d1054b4693

void CollideSweepSinglePrecisionThermalized::run(IBlock *block) {
  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);
  auto force = block->getData<field::GhostLayerField<float, 3>>(forceID);

  auto &omega_even = this->omega_even_;
  auto &seed = this->seed_;
  auto block_offset_0 = this->block_offset_0_;
  auto &omega_odd = this->omega_odd_;
  auto block_offset_1 = this->block_offset_1_;
  auto &kT = this->kT_;
  auto &omega_bulk = this->omega_bulk_;
  auto &time_step = this->time_step_;
  auto &omega_shear = this->omega_shear_;
  auto block_offset_2 = this->block_offset_2_;
  block_offset_generator(block, block_offset_0, block_offset_1, block_offset_2);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force->nrOfGhostLayers()));
  float *RESTRICT const _data_force = force->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()));
  float *RESTRICT _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(), int64_t(cell_idx_c(force->xSize()) + 0));
  const int64_t _size_force_0 = int64_t(cell_idx_c(force->xSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(), int64_t(cell_idx_c(force->ySize()) + 0));
  const int64_t _size_force_1 = int64_t(cell_idx_c(force->ySize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(), int64_t(cell_idx_c(force->zSize()) + 0));
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
  internal_69764eed2d0964e29e3b97d1054b4693::collidesweepsingleprecisionthermalized_collidesweepsingleprecisionthermalized(_data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, block_offset_0, block_offset_1, block_offset_2, kT, omega_bulk, omega_even, omega_odd, omega_shear, seed, time_step);
}

void CollideSweepSinglePrecisionThermalized::runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks, const CellInterval &globalCellInterval, cell_idx_t ghostLayers, IBlock *block) {
  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);
  auto force = block->getData<field::GhostLayerField<float, 3>>(forceID);

  auto &omega_even = this->omega_even_;
  auto &seed = this->seed_;
  auto block_offset_0 = this->block_offset_0_;
  auto &omega_odd = this->omega_odd_;
  auto block_offset_1 = this->block_offset_1_;
  auto &kT = this->kT_;
  auto &omega_bulk = this->omega_bulk_;
  auto &time_step = this->time_step_;
  auto &omega_shear = this->omega_shear_;
  auto block_offset_2 = this->block_offset_2_;
  block_offset_generator(block, block_offset_0, block_offset_1, block_offset_2);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(force->nrOfGhostLayers()));
  float *RESTRICT const _data_force = force->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
  float *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->xSizeWithGhostLayer(), int64_t(cell_idx_c(ci.xSize()) + 0));
  const int64_t _size_force_0 = int64_t(cell_idx_c(ci.xSize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->ySizeWithGhostLayer(), int64_t(cell_idx_c(ci.ySize()) + 0));
  const int64_t _size_force_1 = int64_t(cell_idx_c(ci.ySize()) + 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(force->zSizeWithGhostLayer(), int64_t(cell_idx_c(ci.zSize()) + 0));
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
  internal_69764eed2d0964e29e3b97d1054b4693::collidesweepsingleprecisionthermalized_collidesweepsingleprecisionthermalized(_data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, block_offset_0, block_offset_1, block_offset_2, kT, omega_bulk, omega_even, omega_odd, omega_shear, seed, time_step);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif