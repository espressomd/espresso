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
//! \\file CollideSweepSinglePrecision.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

#include <cmath>

#include "CollideSweepSinglePrecision.h"
#include "core/DataTypes.h"
#include "core/Macros.h"

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

namespace internal_ac3754ab06497f05a231a858fd02d610 {
static FUNC_PREFIX void collidesweepsingleprecision_collidesweepsingleprecision(
    float *RESTRICT const _data_force, float *RESTRICT _data_pdfs,
    int64_t const _size_force_0, int64_t const _size_force_1,
    int64_t const _size_force_2, int64_t const _stride_force_0,
    int64_t const _stride_force_1, int64_t const _stride_force_2,
    int64_t const _stride_force_3, int64_t const _stride_pdfs_0,
    int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
    int64_t const _stride_pdfs_3, float omega_bulk, float omega_even,
    float omega_odd, float omega_shear) {
  const float xi_29 = omega_bulk * 0.50000000000000000f;
  const float xi_56 = omega_shear * 0.041666666666666667f;
  const float xi_61 = omega_bulk * 0.041666666666666667f;
  const float xi_72 = omega_shear * 0.12500000000000000f;
  const float xi_128 = omega_odd * 0.25000000000000000f;
  const float xi_134 = omega_odd * 0.083333333333333333f;
  const float xi_171 = omega_shear * 0.25000000000000000f;
  const float xi_194 = omega_odd * 0.041666666666666667f;
  const float xi_196 = omega_odd * 0.12500000000000000f;
  const float rr_0 = 0.0f;
  const float xi_54 = rr_0 * 0.041666666666666667f;
  const float xi_140 = rr_0 * 0.16666666666666667f;
  const float xi_176 = rr_0 * 0.083333333333333333f;
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
        const float xi_220 = _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0];
        const float xi_221 = _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0];
        const float xi_222 = _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0];
        const float xi_223 = _data_force_20_30_10[_stride_force_0 * ctr_0];
        const float xi_224 = _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0];
        const float xi_225 = _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0];
        const float xi_226 = _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0];
        const float xi_227 = _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0];
        const float xi_228 = _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0];
        const float xi_229 = _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0];
        const float xi_230 = _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0];
        const float xi_231 = _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0];
        const float xi_232 = _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0];
        const float xi_233 = _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0];
        const float xi_234 = _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0];
        const float xi_235 = _data_force_20_31_10[_stride_force_0 * ctr_0];
        const float xi_236 = _data_force_20_32_10[_stride_force_0 * ctr_0];
        const float xi_237 = _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0];
        const float xi_238 = _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0];
        const float xi_239 = _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0];
        const float xi_240 = _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0];
        const float xi_241 = _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0];
        const float xi_0 = xi_226 + xi_231;
        const float xi_1 = xi_0 + xi_222;
        const float xi_2 = xi_220 + xi_227 + xi_233;
        const float xi_3 = xi_234 + xi_238;
        const float xi_4 = xi_221 + xi_232;
        const float xi_5 = xi_224 + xi_237;
        const float xi_6 = xi_229 + xi_241;
        const float xi_9 = -xi_221;
        const float xi_10 = -xi_225 + xi_9;
        const float xi_11 = -xi_229;
        const float xi_12 = -xi_239;
        const float xi_13 = -xi_232;
        const float xi_14 = xi_11 + xi_12 + xi_13;
        const float xi_15 = -xi_237;
        const float xi_16 = -xi_240;
        const float xi_17 = xi_15 + xi_16;
        const float xi_18 = -xi_224;
        const float xi_19 = -xi_234;
        const float xi_20 = xi_18 + xi_19;
        const float xi_21 = -xi_231;
        const float xi_22 = xi_11 + xi_21;
        const float xi_23 = -xi_233;
        const float xi_24 = -xi_241;
        const float xi_25 = xi_18 + xi_227 + xi_23 + xi_24;
        const float xi_30 = xi_235 * 0.16666666666666667f;
        const float xi_31 = xi_235 * 0.083333333333333333f;
        const float xi_43 = xi_223 * 0.16666666666666667f;
        const float xi_44 = xi_223 * 0.083333333333333333f;
        const float xi_50 = xi_236 * 0.16666666666666667f;
        const float xi_51 = xi_236 * 0.083333333333333333f;
        const float xi_68 = xi_235 * 0.25000000000000000f;
        const float xi_73 = xi_235 * xi_72;
        const float xi_106 = -xi_230;
        const float xi_107 = xi_106 + xi_238 * 3.0f + xi_241 * 3.0f;
        const float xi_108 =
            omega_even *
            (xi_107 + xi_220 * 3.0f + xi_224 * -3.0f + xi_227 * -3.0f +
             xi_233 * -3.0f + xi_234 * -3.0f + xi_237 * 3.0f);
        const float xi_109 =
            xi_224 * 2.0f + xi_227 * 2.0f + xi_233 * 2.0f + xi_234 * 2.0f;
        const float xi_110 = xi_109 + xi_222 * 5.0f + xi_232 * 5.0f;
        const float xi_111 =
            omega_even *
            (xi_107 + xi_110 + xi_220 * -2.0f + xi_226 * -5.0f +
             xi_229 * -5.0f + xi_231 * -5.0f + xi_237 * -2.0f + xi_239 * -5.0f);
        const float xi_114 = -xi_227;
        const float xi_115 = xi_114 + xi_19;
        const float xi_116 = -xi_228;
        const float xi_119 = -xi_226;
        const float xi_120 = xi_119 + xi_12 + xi_16 + xi_22;
        const float xi_122 = xi_239 * 2.0f;
        const float xi_123 = xi_226 * 2.0f;
        const float xi_124 = xi_229 * 2.0f + xi_231 * 2.0f;
        const float xi_125 =
            omega_even *
            (xi_106 + xi_110 + xi_122 + xi_123 + xi_124 + xi_220 * 5.0f +
             xi_221 * -7.0f + xi_225 * -7.0f + xi_228 * -7.0f + xi_237 * 5.0f +
             xi_238 * -4.0f + xi_240 * -7.0f + xi_241 * -4.0f);
        const float xi_126 = xi_114 + xi_234;
        const float xi_127 = xi_126 + xi_15 + xi_220 + xi_224 + xi_23;
        const float xi_129 = xi_127 * xi_128;
        const float xi_130 = xi_225 * 2.0f;
        const float xi_131 = xi_240 * 2.0f;
        const float xi_132 = xi_221 * 2.0f + xi_228 * -2.0f;
        const float xi_133 = -xi_130 + xi_131 + xi_132 + xi_15 + xi_2 + xi_20;
        const float xi_135 = xi_133 * xi_134;
        const float xi_136 = -xi_135;
        const float xi_138 = xi_116 + xi_240;
        const float xi_142 = xi_229 + xi_239;
        const float xi_146 = xi_125 * -0.019841269841269841f;
        const float xi_154 = xi_119 + xi_239;
        const float xi_155 = xi_13 + xi_154 + xi_21 + xi_222 + xi_229;
        const float xi_156 = xi_128 * xi_155;
        const float xi_157 = xi_1 + xi_130 - xi_131 + xi_132 + xi_14;
        const float xi_158 = xi_134 * xi_157;
        const float xi_160 = -xi_158;
        const float xi_161 = xi_224 + xi_233;
        const float xi_162 = xi_115 + xi_161 + xi_238 + xi_24;
        const float xi_163 = xi_128 * xi_162;
        const float xi_166 = -xi_122 - xi_123 + xi_124 + xi_25 + xi_3;
        const float xi_167 = xi_134 * xi_166;
        const float xi_168 = -xi_167;
        const float xi_170 = xi_167;
        const float xi_174 = xi_125 * 0.013888888888888889f;
        const float xi_190 = xi_111 * -0.0071428571428571429f;
        const float xi_192 = xi_108 * 0.025000000000000000f;
        const float xi_195 = xi_166 * xi_194;
        const float xi_197 = xi_162 * xi_196;
        const float xi_198 = xi_125 * -0.0039682539682539683f;
        const float xi_202 = xi_133 * xi_194;
        const float xi_203 = xi_127 * xi_196;
        const float xi_209 = xi_111 * 0.017857142857142857f;
        const float xi_212 = xi_155 * xi_196;
        const float xi_213 = xi_157 * xi_194;
        const float xi_32 = rr_0 * xi_31;
        const float xi_45 = rr_0 * xi_44;
        const float xi_52 = rr_0 * xi_51;
        const float xi_55 = xi_223 * xi_54;
        const float xi_60 = xi_235 * xi_54;
        const float xi_82 = xi_236 * xi_54;
        const float vel0Term = xi_1 + xi_228 + xi_240;
        const float vel1Term = xi_2 + xi_225;
        const float vel2Term = xi_239 + xi_3;
        const float rho =
            vel0Term + vel1Term + vel2Term + xi_230 + xi_4 + xi_5 + xi_6;
        const float xi_7 = 1 / (rho);
        const float xi_8 = xi_7 * 0.50000000000000000f;
        const float u_0 = xi_223 * xi_8 + xi_7 * (vel0Term + xi_10 + xi_14);
        const float xi_26 = u_0 * xi_223;
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
        const float xi_96 = u_0 * xi_236;
        const float xi_97 = xi_96 * 0.25000000000000000f;
        const float xi_100 = xi_72 * xi_96;
        const float xi_112 = rho * (u_0 * u_0);
        const float xi_151 = rho * u_0;
        const float xi_152 = -vel0Term + xi_142 + xi_151 + xi_225 + xi_4;
        const float xi_153 = xi_140 * xi_152;
        const float xi_180 = xi_152 * xi_176;
        const float u_1 =
            xi_235 * xi_8 + xi_7 * (vel1Term + xi_17 + xi_20 + xi_228 + xi_9);
        const float xi_27 = u_1 * xi_235;
        const float xi_33 = xi_27 * 0.16666666666666667f;
        const float xi_46 = xi_27 * 0.083333333333333333f;
        const float xi_47 = omega_shear * xi_46;
        const float xi_48 = -xi_33 + xi_47;
        const float xi_63 = -xi_27 * xi_61;
        const float xi_70 = u_1 * 0.25000000000000000f;
        const float xi_71 = xi_223 * xi_70;
        const float xi_75 = u_1 * xi_72;
        const float xi_76 = xi_223 * xi_75;
        const float xi_77 = -xi_69 - xi_71 + xi_74 + xi_76;
        const float xi_79 = xi_69 + xi_71 - xi_74 - xi_76;
        const float xi_87 = xi_236 * xi_70;
        const float xi_89 = xi_236 * xi_75;
        const float xi_94 = -xi_46;
        const float xi_117 = rho * (u_1 * u_1);
        const float xi_118 = xi_10 + xi_116 + xi_117;
        const float xi_137 = rho * u_1;
        const float xi_139 =
            -vel1Term + xi_137 + xi_138 + xi_221 + xi_234 + xi_5;
        const float xi_141 = xi_139 * xi_140;
        const float xi_172 = xi_171 * (u_0 * xi_137 + xi_138 + xi_225 + xi_9);
        const float xi_177 = xi_139 * xi_176;
        const float xi_178 = xi_177;
        const float xi_179 = xi_135 + xi_178;
        const float xi_188 = -xi_177;
        const float xi_189 = xi_136 + xi_188;
        const float xi_204 = xi_178 - xi_202 + xi_203;
        const float xi_205 = xi_188 + xi_202 - xi_203;
        const float u_2 =
            xi_236 * xi_8 + xi_7 * (vel2Term + xi_22 + xi_226 + xi_25);
        const float xi_28 = u_2 * xi_236;
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
        const float xi_98 = u_2 * xi_223;
        const float xi_99 = xi_98 * 0.25000000000000000f;
        const float xi_101 = xi_72 * xi_98;
        const float xi_102 = xi_100 + xi_101 - xi_97 - xi_99;
        const float xi_103 = -xi_100 - xi_101 + xi_97 + xi_99;
        const float xi_104 = -xi_51 + xi_82 + xi_83;
        const float xi_105 = xi_104 + xi_47 + xi_62 + xi_63 + xi_64 + xi_94;
        const float xi_113 = rho * (u_2 * u_2);
        const float xi_121 = omega_bulk * (xi_112 + xi_113 + xi_115 + xi_118 +
                                           xi_120 + xi_18 + xi_23 + xi_230);
        const float xi_143 = -xi_113 + xi_238 + xi_241;
        const float xi_144 =
            omega_shear * (xi_0 + xi_118 + xi_142 + xi_143 + xi_17 - xi_220);
        const float xi_145 = xi_144 * 0.12500000000000000f;
        const float xi_147 =
            omega_shear *
            (xi_10 + xi_109 + xi_112 * 2.0f + xi_116 - xi_117 + xi_120 +
             xi_143 + xi_220 + xi_222 * -2.0f + xi_232 * -2.0f + xi_237);
        const float xi_148 = xi_147 * -0.041666666666666667f;
        const float xi_149 = xi_108 * -0.050000000000000000f + xi_148;
        const float xi_150 =
            xi_111 * 0.014285714285714286f + xi_145 + xi_146 + xi_149;
        const float xi_159 = xi_111 * -0.035714285714285714f + xi_146 +
                             xi_147 * 0.083333333333333333f;
        const float xi_164 =
            rho * u_2 - vel2Term + xi_114 + xi_119 + xi_161 + xi_231 + xi_6;
        const float xi_165 = xi_140 * xi_164;
        const float xi_169 = xi_111 * -0.021428571428571429f +
                             xi_125 * 0.015873015873015873f - xi_145 + xi_149;
        const float xi_173 = xi_144 * 0.062500000000000000f;
        const float xi_175 = -xi_172 + xi_173 + xi_174;
        const float xi_181 = xi_121 * 0.041666666666666667f;
        const float xi_182 = xi_147 * 0.020833333333333333f + xi_181;
        const float xi_183 = -xi_180 + xi_182;
        const float xi_184 = xi_160 + xi_183;
        const float xi_185 = xi_172 + xi_173 + xi_174;
        const float xi_186 = xi_180 + xi_182;
        const float xi_187 = xi_158 + xi_186;
        const float xi_191 = xi_171 * (u_2 * xi_137 + xi_126 + xi_18 + xi_233);
        const float xi_193 = xi_148 + xi_181 + xi_190 + xi_191 + xi_192;
        const float xi_199 = xi_164 * xi_176;
        const float xi_200 = xi_198 + xi_199;
        const float xi_201 = -xi_195 + xi_197 + xi_200;
        const float xi_206 = xi_148 + xi_181 + xi_190 - xi_191 + xi_192;
        const float xi_207 = xi_171 * (u_2 * xi_151 + xi_11 + xi_154 + xi_231);
        const float xi_208 = -xi_173;
        const float xi_210 = -xi_207 + xi_208 + xi_209;
        const float xi_211 = xi_170 + xi_200;
        const float xi_214 = xi_183 - xi_212 + xi_213;
        const float xi_215 = xi_207 + xi_208 + xi_209;
        const float xi_216 = xi_186 + xi_212 - xi_213;
        const float xi_217 = xi_198 - xi_199;
        const float xi_218 = xi_195 - xi_197 + xi_217;
        const float xi_219 = xi_168 + xi_217;
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
            forceTerm_0 + xi_108 * 0.10000000000000000f +
            xi_111 * 0.042857142857142857f + xi_121 * -0.50000000000000000f +
            xi_125 * 0.023809523809523810f + xi_230;
        _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_1 - xi_129 + xi_136 + xi_141 + xi_150 + xi_220;
        _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_2 + xi_129 + xi_135 - xi_141 + xi_150 + xi_237;
        _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_3 - xi_153 + xi_156 + xi_158 + xi_159 + xi_232;
        _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_4 + xi_153 - xi_156 + xi_159 + xi_160 + xi_222;
        _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_5 - xi_163 + xi_165 + xi_168 + xi_169 + xi_238;
        _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_6 + xi_163 - xi_165 + xi_169 + xi_170 + xi_241;
        _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_7 + xi_175 + xi_179 + xi_184 + xi_225;
        _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_8 + xi_179 + xi_185 + xi_187 + xi_228;
        _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_9 + xi_184 + xi_185 + xi_189 + xi_221;
        _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_10 + xi_175 + xi_187 + xi_189 + xi_240;
        _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_11 + xi_193 + xi_201 + xi_204 + xi_227;
        _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_12 + xi_201 + xi_205 + xi_206 + xi_234;
        _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_13 + xi_210 + xi_211 + xi_214 + xi_239;
        _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_14 + xi_211 + xi_215 + xi_216 + xi_226;
        _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_15 + xi_204 + xi_206 + xi_218 + xi_233;
        _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_16 + xi_193 + xi_205 + xi_218 + xi_224;
        _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_17 + xi_214 + xi_215 + xi_219 + xi_229;
        _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_18 + xi_210 + xi_216 + xi_219 + xi_231;
      }
    }
  }
}
} // namespace internal_ac3754ab06497f05a231a858fd02d610

void CollideSweepSinglePrecision::run(IBlock *block) {
  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);
  auto force = block->getData<field::GhostLayerField<float, 3>>(forceID);

  auto &omega_odd = this->omega_odd_;
  auto &omega_even = this->omega_even_;
  auto &omega_shear = this->omega_shear_;
  auto &omega_bulk = this->omega_bulk_;
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
  internal_ac3754ab06497f05a231a858fd02d610::
      collidesweepsingleprecision_collidesweepsingleprecision(
          _data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2,
          _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3,
          _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3,
          omega_bulk, omega_even, omega_odd, omega_shear);
}

void CollideSweepSinglePrecision::runOnCellInterval(
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
  auto &omega_even = this->omega_even_;
  auto &omega_shear = this->omega_shear_;
  auto &omega_bulk = this->omega_bulk_;
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
  internal_ac3754ab06497f05a231a858fd02d610::
      collidesweepsingleprecision_collidesweepsingleprecision(
          _data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2,
          _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3,
          _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3,
          omega_bulk, omega_even, omega_odd, omega_shear);
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