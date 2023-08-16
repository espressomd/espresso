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
//! \\file CollideSweepDoublePrecisionThermalized.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

// kernel generated with pystencils v1.3.1+2.g60e24c4, lbmpy v1.3.1+6.gcd1bc2f.dirty, lbmpy_walberla/pystencils_walberla from waLBerla commit 065ce5f311850371a97ac4766f47dbb5ca8424ba

#include <cmath>

#include "CollideSweepDoublePrecisionThermalized.h"
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

namespace internal_0d943397135d13b4628c5752888935d7 {
static FUNC_PREFIX void collidesweepdoubleprecisionthermalized_collidesweepdoubleprecisionthermalized(double *RESTRICT const _data_force, double *RESTRICT _data_pdfs, int64_t const _size_force_0, int64_t const _size_force_1, int64_t const _size_force_2, int64_t const _stride_force_0, int64_t const _stride_force_1, int64_t const _stride_force_2, int64_t const _stride_force_3, int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3, uint32_t block_offset_0, uint32_t block_offset_1, uint32_t block_offset_2, double kT, double omega_bulk, double omega_even, double omega_odd, double omega_shear, uint32_t seed, uint32_t time_step) {
  const double xi_21 = omega_bulk * 0.5;
  const double xi_48 = omega_shear * 0.041666666666666664;
  const double xi_52 = omega_bulk * 0.041666666666666664;
  const double xi_63 = omega_shear * 0.125;
  const double xi_98 = 3.7416573867739413;
  const double xi_101 = 5.4772255750516612;
  const double xi_105 = 2.4494897427831779;
  const double xi_108 = 8.3666002653407556;
  const double xi_151 = omega_odd * 0.25;
  const double xi_161 = omega_odd * 0.083333333333333329;
  const double xi_174 = 1.7320508075688772;
  const double xi_218 = omega_shear * 0.25;
  const double xi_224 = omega_odd * 0.041666666666666664;
  const double xi_227 = omega_odd * 0.125;
  const double rr_0 = 0.0;
  const double xi_46 = rr_0 * 0.041666666666666664;
  for (int64_t ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_20_37 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_315 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_314 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_316 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_34 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_31 = _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_32 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    double *RESTRICT _data_pdfs_20_36 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_33 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    double *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    double *RESTRICT _data_pdfs_20_313 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_317 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_39 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    double *RESTRICT _data_force_20_31 = _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    double *RESTRICT _data_pdfs_20_310 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    double *RESTRICT _data_force_20_32 = _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    double *RESTRICT _data_pdfs_20_318 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_35 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_312 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_311 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_38 = _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_20_37_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_37;
      double *RESTRICT _data_pdfs_20_315_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_315;
      double *RESTRICT _data_pdfs_20_314_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_314;
      double *RESTRICT _data_pdfs_20_316_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_316;
      double *RESTRICT _data_pdfs_20_34_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      double *RESTRICT _data_pdfs_20_31_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_31;
      double *RESTRICT _data_pdfs_20_32_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_32;
      double *RESTRICT _data_pdfs_20_30_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      double *RESTRICT _data_pdfs_20_36_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_36;
      double *RESTRICT _data_pdfs_20_33_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      double *RESTRICT _data_force_20_30_10 = _stride_force_1 * ctr_1 + _data_force_20_30;
      double *RESTRICT _data_pdfs_20_313_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_313;
      double *RESTRICT _data_pdfs_20_317_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_317;
      double *RESTRICT _data_pdfs_20_39_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_39;
      double *RESTRICT _data_force_20_31_10 = _stride_force_1 * ctr_1 + _data_force_20_31;
      double *RESTRICT _data_pdfs_20_310_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_310;
      double *RESTRICT _data_force_20_32_10 = _stride_force_1 * ctr_1 + _data_force_20_32;
      double *RESTRICT _data_pdfs_20_318_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_318;
      double *RESTRICT _data_pdfs_20_35_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_35;
      double *RESTRICT _data_pdfs_20_312_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_312;
      double *RESTRICT _data_pdfs_20_311_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_311;
      double *RESTRICT _data_pdfs_20_38_10 = _stride_pdfs_1 * ctr_1 + _data_pdfs_20_38;
      for (int64_t ctr_0 = 0; ctr_0 < _size_force_0; ctr_0 += 1) {
        const double xi_264 = _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0];
        const double xi_265 = _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0];
        const double xi_266 = _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0];
        const double xi_267 = _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0];
        const double xi_268 = _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0];
        const double xi_269 = _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0];
        const double xi_270 = _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0];
        const double xi_271 = _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0];
        const double xi_272 = _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0];
        const double xi_273 = _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0];
        const double xi_274 = _data_force_20_30_10[_stride_force_0 * ctr_0];
        const double xi_275 = _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0];
        const double xi_276 = _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0];
        const double xi_277 = _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0];
        const double xi_278 = _data_force_20_31_10[_stride_force_0 * ctr_0];
        const double xi_279 = _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0];
        const double xi_280 = _data_force_20_32_10[_stride_force_0 * ctr_0];
        const double xi_281 = _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0];
        const double xi_282 = _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0];
        const double xi_283 = _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0];
        const double xi_284 = _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0];
        const double xi_285 = _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0];

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
        const double xi_2 = xi_279 + xi_285;
        const double xi_3 = xi_2 + xi_266 + xi_281;
        const double xi_4 = xi_264 + xi_269;
        const double xi_5 = xi_265 + xi_284;
        const double xi_6 = xi_282 + xi_283;
        const double xi_7 = xi_273 + xi_276;
        const double xi_8 = xi_277 + xi_7;
        const double xi_9 = xi_267 + xi_270;
        const double xi_12 = xi_264 + xi_275;
        const double xi_13 = xi_279 + xi_285 * -1.0;
        const double xi_14 = xi_277 + xi_283;
        const double xi_15 = xi_284 * -1.0;
        const double xi_16 = xi_266 * -1.0;
        const double xi_17 = xi_267 + xi_276;
        const double xi_22 = xi_278 * 0.16666666666666666;
        const double xi_23 = xi_278 * 0.083333333333333329;
        const double xi_34 = xi_274 * 0.16666666666666666;
        const double xi_35 = xi_274 * 0.083333333333333329;
        const double xi_40 = xi_280 * 0.16666666666666666;
        const double xi_41 = xi_280 * 0.083333333333333329;
        const double xi_59 = xi_278 * 0.25;
        const double xi_64 = xi_278 * xi_63;
        const double xi_99 = random_6_0 - 0.5;
        const double xi_102 = random_7_0 - 0.5;
        const double xi_104 = random_2_1 - 0.5;
        const double xi_109 = random_6_1 - 0.5;
        const double xi_113 = xi_271 * -1.0;
        const double xi_119 = xi_12 + xi_3;
        const double xi_122 = xi_265 * 2.0 + xi_267 * 2.0 + xi_268 * 5.0 + xi_273 * 5.0 + xi_283 * 2.0 + xi_284 * 2.0;
        const double xi_125 = xi_275 * 2.0;
        const double xi_126 = xi_266 * 2.0;
        const double xi_127 = xi_276 * 2.0 + xi_281 * 2.0;
        const double xi_130 = xi_269 * -1.0;
        const double xi_134 = random_0_1 - 0.5;
        const double xi_138 = xi_276 * -1.0;
        const double xi_139 = xi_281 * -1.0;
        const double xi_140 = xi_282 * -1.0;
        const double xi_141 = xi_272 * -1.0;
        const double xi_146 = xi_265 * -1.0;
        const double xi_147 = xi_15 + xi_283;
        const double xi_148 = xi_270 * -1.0;
        const double xi_149 = xi_148 + xi_269;
        const double xi_150 = xi_146 + xi_147 + xi_149 + xi_267;
        const double xi_152 = random_5_1 - 0.5;
        const double xi_156 = xi_267 * -1.0;
        const double xi_157 = xi_283 * -1.0;
        const double xi_158 = xi_285 * 2.0;
        const double xi_159 = xi_264 * -2.0 + xi_279 * 2.0;
        const double xi_160 = xi_149 + xi_156 + xi_157 + xi_158 * -1.0 + xi_159 + xi_277 * 2.0 + xi_5;
        const double xi_162 = xi_160 * xi_161;
        const double xi_163 = random_3_0 - 0.5;
        const double xi_175 = random_0_0 - 0.5;
        const double xi_186 = xi_16 + xi_275;
        const double xi_187 = xi_139 + xi_186;
        const double xi_188 = xi_187 + xi_268 + xi_273 * -1.0 + xi_276;
        const double xi_189 = random_4_1 - 0.5;
        const double xi_190 = xi_268 * -1.0;
        const double xi_191 = xi_158 * -1.0 + xi_159 * -1.0 + xi_187 * -1.0 + xi_190 * -1.0 + xi_277 * 2.0 + xi_7 * -1.0;
        const double xi_192 = xi_161 * xi_191;
        const double xi_193 = random_4_0 - 0.5;
        const double xi_199 = xi_146 + xi_156 + xi_284;
        const double xi_200 = xi_140 * -1.0 + xi_199 * -1.0 + xi_272 * -1.0 + xi_283 * -1.0;
        const double xi_201 = random_5_0 - 0.5;
        const double xi_202 = xi_125 * -1.0 + xi_126 * -1.0 + xi_127 + xi_141 + xi_199 + xi_6;
        const double xi_203 = xi_161 * xi_202;
        const double xi_204 = random_3_1 - 0.5;
        const double xi_225 = xi_160 * xi_224;
        const double xi_228 = xi_150 * xi_227;
        const double xi_234 = xi_202 * xi_224;
        const double xi_235 = xi_200 * xi_227;
        const double xi_255 = xi_188 * xi_227;
        const double xi_256 = xi_191 * xi_224;
        const double xi_24 = rr_0 * xi_23;
        const double xi_36 = rr_0 * xi_35;
        const double xi_42 = rr_0 * xi_41;
        const double xi_47 = xi_274 * xi_46;
        const double xi_51 = xi_278 * xi_46;
        const double xi_73 = xi_280 * xi_46;
        const double vel0Term = xi_268 + xi_3;
        const double vel1Term = xi_4 + xi_5;
        const double vel2Term = xi_275 + xi_6;
        const double delta_rho = vel0Term + vel1Term + vel2Term + xi_271 + xi_272 + xi_8 + xi_9;
        const double xi_10 = delta_rho + 1.0;
        const double xi_96 = kT * xi_10;
        const double xi_97 = pow(xi_96 * (-1.0 * ((omega_even * -1.0 + 1.0) * (omega_even * -1.0 + 1.0)) + 1.0), 0.5);
        const double xi_100 = xi_97 * xi_98 * xi_99;
        const double xi_103 = xi_101 * xi_102 * xi_97;
        const double xi_106 = pow(xi_96 * (-1.0 * ((omega_bulk * -1.0 + 1.0) * (omega_bulk * -1.0 + 1.0)) + 1.0), 0.5);
        const double xi_107 = xi_104 * xi_105 * xi_106;
        const double xi_110 = xi_108 * xi_109 * xi_97;
        const double xi_132 = xi_100 * 0.11904761904761904;
        const double xi_135 = pow(xi_96 * (-1.0 * ((omega_shear * -1.0 + 1.0) * (omega_shear * -1.0 + 1.0)) + 1.0), 0.5);
        const double xi_136 = xi_135 * 0.5;
        const double xi_137 = xi_134 * xi_136;
        const double xi_153 = pow(xi_96 * (-1.0 * ((omega_odd * -1.0 + 1.0) * (omega_odd * -1.0 + 1.0)) + 1.0), 0.5);
        const double xi_154 = xi_153 * 1.4142135623730951;
        const double xi_155 = xi_154 * 0.5;
        const double xi_164 = xi_105 * xi_153;
        const double xi_165 = xi_164 * 0.16666666666666666;
        const double xi_166 = xi_163 * xi_165;
        const double xi_167 = xi_162 + xi_166;
        const double xi_168 = xi_150 * xi_151 + xi_152 * xi_155 + xi_167;
        const double xi_170 = xi_103 * 0.10000000000000001;
        const double xi_176 = xi_135 * xi_174 * xi_175;
        const double xi_177 = xi_176 * 0.16666666666666666;
        const double xi_185 = xi_110 * 0.071428571428571425;
        const double xi_194 = xi_165 * xi_193;
        const double xi_195 = xi_192 + xi_194;
        const double xi_196 = xi_151 * xi_188 + xi_155 * xi_189 + xi_195;
        const double xi_198 = xi_110 * 0.042857142857142858;
        const double xi_205 = xi_165 * xi_204;
        const double xi_206 = xi_203 + xi_205;
        const double xi_207 = xi_151 * xi_200 + xi_155 * xi_201 + xi_206;
        const double xi_208 = xi_134 * xi_135 * 0.25;
        const double xi_211 = xi_100 * 0.083333333333333329;
        const double xi_215 = xi_192 * -1.0 + xi_194 * -1.0;
        const double xi_216 = xi_136 * (random_1_0 - 0.5);
        const double xi_223 = xi_136 * (random_2_0 - 0.5);
        const double xi_229 = xi_164 * 0.083333333333333329;
        const double xi_230 = xi_163 * xi_229;
        const double xi_231 = xi_154 * 0.25;
        const double xi_232 = xi_152 * xi_231;
        const double xi_236 = xi_204 * xi_229;
        const double xi_237 = xi_201 * xi_231;
        const double xi_238 = xi_234 * -1.0 + xi_235 + xi_236 * -1.0 + xi_237;
        const double xi_240 = xi_110 * 0.014285714285714285;
        const double xi_242 = xi_100 * 0.023809523809523808;
        const double xi_245 = xi_234 + xi_235 * -1.0 + xi_236 + xi_237 * -1.0;
        const double xi_247 = xi_208 * -1.0;
        const double xi_250 = xi_110 * 0.035714285714285712;
        const double xi_252 = xi_136 * (random_1_1 - 0.5);
        const double xi_257 = xi_189 * xi_231;
        const double xi_258 = xi_193 * xi_229;
        const double xi_259 = xi_255 * -1.0 + xi_256 + xi_257 * -1.0 + xi_258;
        const double xi_261 = xi_255 + xi_256 * -1.0 + xi_257 + xi_258 * -1.0;
        const double rho = xi_10;
        const double xi_0 = ((1.0) / (rho));
        const double xi_11 = xi_0 * 0.5;
        const double u_0 = xi_0 * (vel0Term + xi_12 * -1.0 + xi_8 * -1.0) + xi_11 * xi_274;
        const double xi_18 = u_0 * xi_274;
        const double xi_29 = xi_18 * 0.16666666666666666;
        const double xi_30 = xi_29 * -1.0;
        const double xi_31 = xi_18 * 0.083333333333333329;
        const double xi_32 = omega_shear * xi_31 + xi_30;
        const double xi_49 = xi_18 * xi_48 + xi_30;
        const double xi_50 = xi_35 + xi_47 * -1.0 + xi_49;
        const double xi_53 = xi_18 * xi_52;
        const double xi_60 = u_0 * xi_59;
        const double xi_65 = u_0 * xi_64;
        const double xi_69 = xi_35 * -1.0 + xi_47 + xi_49;
        const double xi_76 = omega_shear * u_0 * xi_274 * -0.083333333333333329;
        const double xi_86 = u_0 * xi_280;
        const double xi_87 = xi_86 * 0.25;
        const double xi_90 = xi_63 * xi_86;
        const double xi_112 = u_0 * u_0;
        const double u_1 = xi_0 * (vel1Term + xi_13 * -1.0 + xi_14 * -1.0 + xi_9 * -1.0) + xi_11 * xi_278;
        const double xi_19 = u_1 * xi_278;
        const double xi_27 = xi_19 * 0.16666666666666666;
        const double xi_37 = omega_shear * u_1 * xi_278 * -0.083333333333333329;
        const double xi_43 = xi_27 * -1.0;
        const double xi_44 = xi_19 * 0.083333333333333329;
        const double xi_54 = xi_19 * xi_52;
        const double xi_61 = u_1 * 0.25;
        const double xi_62 = xi_274 * xi_61;
        const double xi_66 = u_1 * xi_63;
        const double xi_67 = xi_274 * xi_66;
        const double xi_68 = xi_60 + xi_62 + xi_65 * -1.0 + xi_67 * -1.0;
        const double xi_70 = xi_60 * -1.0 + xi_62 * -1.0 + xi_65 + xi_67;
        const double xi_78 = xi_280 * xi_61;
        const double xi_80 = xi_280 * xi_66;
        const double xi_111 = rho * (u_1 * u_1);
        const double xi_118 = xi_111 * -1.0;
        const double xi_217 = rho * u_1;
        const double xi_219 = xi_218 * (u_0 * xi_217 + xi_13 + xi_264 + xi_277 * -1.0);
        const double xi_220 = xi_216 * -1.0 + xi_219 * -1.0;
        const double xi_221 = xi_216 + xi_219;
        const double u_2 = xi_0 * (vel2Term + xi_15 * -1.0 + xi_16 * -1.0 + xi_17 * -1.0 + xi_265 * -1.0 + xi_272 * -1.0 + xi_281 * -1.0) + xi_11 * xi_280;
        const double xi_20 = u_2 * xi_280;
        const double xi_25 = xi_20 * 0.16666666666666666;
        const double xi_26 = xi_25 * -1.0;
        const double xi_28 = xi_20 * 0.083333333333333329;
        const double xi_33 = omega_shear * xi_27 * -1.0 + omega_shear * xi_28 + xi_19 * 0.33333333333333331 + xi_26 + xi_32;
        const double xi_38 = omega_shear * u_2 * xi_280 * -0.083333333333333329;
        const double xi_39 = omega_shear * xi_29 + u_0 * xi_274 * -0.33333333333333331 + xi_25 + xi_27 + xi_37 + xi_38;
        const double xi_45 = omega_shear * xi_25 * -1.0 + omega_shear * xi_44 + xi_20 * 0.33333333333333331 + xi_32 + xi_43;
        const double xi_55 = xi_20 * xi_52;
        const double xi_56 = xi_19 * xi_48 + xi_43 + xi_53 + xi_54 + xi_55;
        const double xi_57 = xi_23 * -1.0 + xi_51 + xi_56;
        const double xi_58 = xi_28 + xi_38 + xi_57;
        const double xi_71 = xi_23 + xi_51 * -1.0 + xi_56;
        const double xi_72 = xi_28 + xi_38 + xi_71;
        const double xi_74 = xi_20 * xi_48 + xi_26;
        const double xi_75 = xi_41 * -1.0 + xi_73 + xi_74;
        const double xi_77 = xi_31 + xi_57 + xi_76;
        const double xi_79 = u_2 * xi_59;
        const double xi_81 = u_2 * xi_64;
        const double xi_82 = xi_78 * -1.0 + xi_79 * -1.0 + xi_80 + xi_81;
        const double xi_83 = xi_31 + xi_71 + xi_76;
        const double xi_84 = xi_78 + xi_79 + xi_80 * -1.0 + xi_81 * -1.0;
        const double xi_85 = xi_37 + xi_44 + xi_53 + xi_54 + xi_55 + xi_75;
        const double xi_88 = u_2 * xi_274;
        const double xi_89 = xi_88 * 0.25;
        const double xi_91 = xi_63 * xi_88;
        const double xi_92 = xi_87 + xi_89 + xi_90 * -1.0 + xi_91 * -1.0;
        const double xi_93 = xi_87 * -1.0 + xi_89 * -1.0 + xi_90 + xi_91;
        const double xi_94 = xi_41 + xi_73 * -1.0 + xi_74;
        const double xi_95 = xi_37 + xi_44 + xi_53 + xi_54 + xi_55 + xi_94;
        const double xi_114 = rho * (u_2 * u_2);
        const double xi_115 = xi_113 + xi_114 * 0.66666666666666663 + xi_272 * 3.0 + xi_282 * 3.0;
        const double xi_116 = rho * xi_112 * 1.6666666666666667 + xi_111 * 0.66666666666666663 + xi_115 + xi_265 * -3.0 + xi_267 * -3.0 + xi_269 * 3.0 + xi_270 * 3.0 + xi_283 * -3.0 + xi_284 * -3.0;
        const double xi_117 = omega_even * xi_116;
        const double xi_120 = rho * xi_112 + xi_113 * -1.0 + xi_114 + xi_118 * -1.0 + xi_119 * -1.0 + xi_14 * -1.0 + xi_17 * -1.0 + xi_5 * -1.0;
        const double xi_121 = omega_bulk * xi_120;
        const double xi_123 = xi_111 * 2.3333333333333335 + xi_115 + xi_122 + xi_266 * -5.0 + xi_269 * -2.0 + xi_270 * -2.0 + xi_275 * -5.0 + xi_276 * -5.0 + xi_281 * -5.0;
        const double xi_124 = omega_even * xi_123;
        const double xi_128 = xi_113 + xi_114 * 3.0 + xi_122 + xi_125 + xi_126 + xi_127 + xi_264 * -7.0 + xi_269 * 5.0 + xi_270 * 5.0 + xi_272 * -4.0 + xi_277 * -7.0 + xi_279 * -7.0 + xi_282 * -4.0 + xi_285 * -7.0;
        const double xi_129 = omega_even * xi_128;
        const double xi_131 = xi_129 * 0.01984126984126984;
        const double xi_133 = xi_131 + xi_132;
        const double xi_142 = xi_114 + xi_140 + xi_141 + xi_277;
        const double xi_143 = omega_shear * (xi_118 * -1.0 + xi_138 * -1.0 + xi_139 * -1.0 + xi_142 * -1.0 + xi_16 * -1.0 + xi_2 * -1.0 + xi_270 * -1.0 + xi_275 + xi_4 * -1.0);
        const double xi_144 = xi_143 * 0.125;
        const double xi_145 = xi_137 * -1.0 + xi_144 * -1.0;
        const double xi_169 = xi_117 * 0.050000000000000003;
        const double xi_171 = rho * xi_112 * 2.0 + xi_111 * -1.0 + xi_119 * -1.0 + xi_130 * -1.0 + xi_142 * -1.0 + xi_148 * -1.0 + xi_265 * 2.0 + xi_267 * 2.0 + xi_268 * -2.0 + xi_273 * -2.0 + xi_276 * -1.0 + xi_283 * 2.0 + xi_284 * 2.0;
        const double xi_172 = omega_shear * xi_171;
        const double xi_173 = xi_172 * 0.041666666666666664;
        const double xi_178 = xi_173 + xi_177;
        const double xi_179 = xi_169 + xi_170 + xi_178;
        const double xi_180 = xi_131 * -1.0 + xi_132 * -1.0;
        const double xi_181 = xi_137 + xi_144;
        const double xi_182 = xi_173 * -1.0 + xi_177 * -1.0;
        const double xi_183 = xi_169 * -1.0 + xi_170 * -1.0 + xi_182;
        const double xi_184 = xi_124 * 0.035714285714285712;
        const double xi_197 = xi_124 * 0.021428571428571429;
        const double xi_209 = xi_143 * 0.0625;
        const double xi_210 = xi_129 * 0.013888888888888888;
        const double xi_212 = xi_107 * 0.083333333333333329 + xi_121 * 0.041666666666666664;
        const double xi_213 = xi_172 * 0.020833333333333332 + xi_176 * 0.083333333333333329 + xi_212;
        const double xi_214 = xi_167 + xi_208 + xi_209 + xi_210 + xi_211 + xi_213;
        const double xi_222 = xi_162 * -1.0 + xi_166 * -1.0 + xi_208 + xi_209 + xi_210 + xi_211 + xi_213;
        const double xi_226 = xi_218 * (u_2 * xi_217 + xi_147 + xi_156 + xi_265);
        const double xi_233 = xi_223 + xi_225 * -1.0 + xi_226 + xi_228 + xi_230 * -1.0 + xi_232;
        const double xi_239 = xi_124 * 0.0071428571428571426;
        const double xi_241 = xi_129 * 0.003968253968253968;
        const double xi_243 = xi_241 * -1.0 + xi_242 * -1.0;
        const double xi_244 = xi_103 * 0.050000000000000003 + xi_117 * 0.025000000000000001 + xi_182 + xi_212 + xi_239 * -1.0 + xi_240 * -1.0 + xi_243;
        const double xi_246 = omega_bulk * xi_120 * -0.041666666666666664 + omega_even * xi_116 * -0.025000000000000001 + xi_101 * xi_102 * xi_97 * -0.050000000000000003 + xi_104 * xi_105 * xi_106 * -0.083333333333333329 + xi_178 + xi_239 + xi_240 + xi_241 + xi_242;
        const double xi_248 = xi_209 * -1.0;
        const double xi_249 = xi_124 * 0.017857142857142856;
        const double xi_251 = xi_206 + xi_213 + xi_243 + xi_247 + xi_248 + xi_249 + xi_250;
        const double xi_253 = xi_218 * (rho * u_0 * u_2 + xi_138 + xi_186 + xi_281);
        const double xi_254 = xi_252 * -1.0 + xi_253 * -1.0;
        const double xi_260 = xi_252 + xi_253;
        const double xi_262 = xi_223 + xi_225 + xi_226 + xi_228 * -1.0 + xi_230 + xi_232 * -1.0;
        const double xi_263 = xi_203 * -1.0 + xi_205 * -1.0 + xi_213 + xi_243 + xi_247 + xi_248 + xi_249 + xi_250;
        const double forceTerm_0 = xi_18 * xi_21 + xi_18 * -1.0 + xi_19 * xi_21 + xi_19 * -1.0 + xi_20 * xi_21 + xi_20 * -1.0;
        const double forceTerm_1 = xi_22 + xi_24 * -1.0 + xi_33;
        const double forceTerm_2 = xi_22 * -1.0 + xi_24 + xi_33;
        const double forceTerm_3 = xi_34 * -1.0 + xi_36 + xi_39 * -1.0;
        const double forceTerm_4 = xi_34 + xi_36 * -1.0 + xi_39 * -1.0;
        const double forceTerm_5 = xi_40 + xi_42 * -1.0 + xi_45;
        const double forceTerm_6 = xi_40 * -1.0 + xi_42 + xi_45;
        const double forceTerm_7 = xi_50 * -1.0 + xi_58 * -1.0 + xi_68 * -1.0;
        const double forceTerm_8 = xi_58 * -1.0 + xi_69 * -1.0 + xi_70 * -1.0;
        const double forceTerm_9 = xi_50 * -1.0 + xi_70 * -1.0 + xi_72 * -1.0;
        const double forceTerm_10 = xi_68 * -1.0 + xi_69 * -1.0 + xi_72 * -1.0;
        const double forceTerm_11 = xi_75 * -1.0 + xi_77 * -1.0 + xi_82 * -1.0;
        const double forceTerm_12 = xi_75 * -1.0 + xi_83 * -1.0 + xi_84 * -1.0;
        const double forceTerm_13 = xi_50 * -1.0 + xi_85 * -1.0 + xi_92 * -1.0;
        const double forceTerm_14 = xi_69 * -1.0 + xi_85 * -1.0 + xi_93 * -1.0;
        const double forceTerm_15 = xi_77 * -1.0 + xi_84 * -1.0 + xi_94 * -1.0;
        const double forceTerm_16 = xi_82 * -1.0 + xi_83 * -1.0 + xi_94 * -1.0;
        const double forceTerm_17 = xi_50 * -1.0 + xi_93 * -1.0 + xi_95 * -1.0;
        const double forceTerm_18 = xi_69 * -1.0 + xi_92 * -1.0 + xi_95 * -1.0;
        _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0] = forceTerm_0 + xi_100 * 0.14285714285714285 + xi_103 * 0.20000000000000001 + xi_107 * -1.0 + xi_110 * 0.085714285714285715 + xi_117 * 0.10000000000000001 + xi_121 * -0.5 + xi_124 * 0.042857142857142858 + xi_129 * 0.023809523809523808 + xi_271;
        _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0] = forceTerm_1 + omega_even * xi_123 * 0.014285714285714285 + xi_108 * xi_109 * xi_97 * 0.028571428571428571 + xi_130 * -1.0 + xi_133 * -1.0 + xi_145 * -1.0 + xi_168 * -1.0 + xi_179 * -1.0;
        _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0] = forceTerm_2 + xi_110 * 0.028571428571428571 + xi_124 * 0.014285714285714285 + xi_168 + xi_180 + xi_181 + xi_183 + xi_270;
        _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0] = forceTerm_3 + xi_172 * 0.083333333333333329 + xi_176 * 0.33333333333333331 + xi_180 + xi_184 * -1.0 + xi_185 * -1.0 + xi_196 + xi_273;
        _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0] = forceTerm_4 + omega_shear * xi_171 * 0.083333333333333329 + xi_133 * -1.0 + xi_135 * xi_174 * xi_175 * 0.33333333333333331 + xi_184 * -1.0 + xi_185 * -1.0 + xi_190 * -1.0 + xi_196 * -1.0;
        _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0] = forceTerm_5 + omega_even * xi_128 * 0.015873015873015872 + xi_140 * -1.0 + xi_179 * -1.0 + xi_181 * -1.0 + xi_197 * -1.0 + xi_198 * -1.0 + xi_207 * -1.0 + xi_97 * xi_98 * xi_99 * 0.095238095238095233;
        _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0] = forceTerm_6 + xi_100 * 0.095238095238095233 + xi_129 * 0.015873015873015872 + xi_145 + xi_183 + xi_197 * -1.0 + xi_198 * -1.0 + xi_207 + xi_272;
        _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0] = forceTerm_7 + xi_214 + xi_215 + xi_220 + xi_264;
        _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0] = forceTerm_8 + xi_195 + xi_214 + xi_221 + xi_285;
        _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0] = forceTerm_9 + xi_215 + xi_221 + xi_222 + xi_277;
        _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0] = forceTerm_10 + xi_195 + xi_220 + xi_222 + xi_279;
        _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0] = forceTerm_11 + xi_233 + xi_238 + xi_244 + xi_284;
        _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0] = forceTerm_12 + xi_157 * -1.0 + xi_233 * -1.0 + xi_245 * -1.0 + xi_246 * -1.0;
        _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0] = forceTerm_13 + xi_251 + xi_254 + xi_259 + xi_275;
        _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0] = forceTerm_14 + xi_251 + xi_260 + xi_261 + xi_266;
        _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0] = forceTerm_15 + xi_146 * -1.0 + xi_238 * -1.0 + xi_246 * -1.0 + xi_262 * -1.0;
        _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0] = forceTerm_16 + xi_244 + xi_245 + xi_262 + xi_267;
        _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0] = forceTerm_17 + xi_259 + xi_260 + xi_263 + xi_276;
        _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0] = forceTerm_18 + xi_254 + xi_261 + xi_263 + xi_281;
      }
    }
  }
}
} // namespace internal_0d943397135d13b4628c5752888935d7

void CollideSweepDoublePrecisionThermalized::run(IBlock *block) {
  auto force = block->getData<field::GhostLayerField<double, 3>>(forceID);
  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);

  auto &omega_odd = this->omega_odd_;
  auto &time_step = this->time_step_;
  auto &omega_shear = this->omega_shear_;
  auto &kT = this->kT_;
  auto block_offset_1 = this->block_offset_1_;
  auto &omega_even = this->omega_even_;
  auto block_offset_2 = this->block_offset_2_;
  auto block_offset_0 = this->block_offset_0_;
  auto &seed = this->seed_;
  auto &omega_bulk = this->omega_bulk_;
  block_offset_generator(block, block_offset_0, block_offset_1, block_offset_2);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force->nrOfGhostLayers()));
  double *RESTRICT const _data_force = force->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()));
  double *RESTRICT _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
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
  internal_0d943397135d13b4628c5752888935d7::collidesweepdoubleprecisionthermalized_collidesweepdoubleprecisionthermalized(_data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, block_offset_0, block_offset_1, block_offset_2, kT, omega_bulk, omega_even, omega_odd, omega_shear, seed, time_step);
}

void CollideSweepDoublePrecisionThermalized::runOnCellInterval(const shared_ptr<StructuredBlockStorage> &blocks, const CellInterval &globalCellInterval, cell_idx_t ghostLayers, IBlock *block) {
  CellInterval ci = globalCellInterval;
  CellInterval blockBB = blocks->getBlockCellBB(*block);
  blockBB.expand(ghostLayers);
  ci.intersect(blockBB);
  blocks->transformGlobalToBlockLocalCellInterval(ci, *block);
  if (ci.empty())
    return;

  auto force = block->getData<field::GhostLayerField<double, 3>>(forceID);
  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);

  auto &omega_odd = this->omega_odd_;
  auto &time_step = this->time_step_;
  auto &omega_shear = this->omega_shear_;
  auto &kT = this->kT_;
  auto block_offset_1 = this->block_offset_1_;
  auto &omega_even = this->omega_even_;
  auto block_offset_2 = this->block_offset_2_;
  auto block_offset_0 = this->block_offset_0_;
  auto &seed = this->seed_;
  auto &omega_bulk = this->omega_bulk_;
  block_offset_generator(block, block_offset_0, block_offset_1, block_offset_2);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(force->nrOfGhostLayers()));
  double *RESTRICT const _data_force = force->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
  double *RESTRICT _data_pdfs = pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
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
  internal_0d943397135d13b4628c5752888935d7::collidesweepdoubleprecisionthermalized_collidesweepdoubleprecisionthermalized(_data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3, block_offset_0, block_offset_1, block_offset_2, kT, omega_bulk, omega_even, omega_odd, omega_shear, seed, time_step);
}

} // namespace pystencils
} // namespace walberla

#if (defined WALBERLA_CXX_COMPILER_IS_GNU) || (defined WALBERLA_CXX_COMPILER_IS_CLANG)
#pragma GCC diagnostic pop
#endif

#if (defined WALBERLA_CXX_COMPILER_IS_INTEL)
#pragma warning pop
#endif