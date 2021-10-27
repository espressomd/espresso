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
//! \\file CollideSweep.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

#include <cmath>

#include "CollideSweep.h"
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

namespace internal_collidesweep {
static FUNC_PREFIX void
collidesweep(double *RESTRICT const _data_force, double *RESTRICT _data_pdfs,
             int64_t const _size_force_0, int64_t const _size_force_1,
             int64_t const _size_force_2, int64_t const _stride_force_0,
             int64_t const _stride_force_1, int64_t const _stride_force_2,
             int64_t const _stride_force_3, int64_t const _stride_pdfs_0,
             int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
             int64_t const _stride_pdfs_3, double omega_bulk, double omega_even,
             double omega_odd, double omega_shear) {
  const double xi_35 = -omega_shear + 2.0;
  const double xi_36 = xi_35 * 0.5;
  const double xi_41 = xi_35 * 0.0833333333333333;
  const double xi_46 = xi_35 * 0.166666666666667;
  const double xi_56 = xi_35 * 0.25;
  const double xi_61 = xi_35 * 0.0416666666666667;
  const double xi_106 = omega_odd * 0.25;
  const double xi_112 = omega_odd * 0.0833333333333333;
  const double xi_149 = omega_shear * 0.25;
  const double xi_172 = omega_odd * 0.0416666666666667;
  const double xi_174 = omega_odd * 0.125;
  const int64_t rr_0 = 0.0;
  const double xi_118 = rr_0 * 0.166666666666667;
  const double xi_154 = rr_0 * 0.0833333333333333;
  for (int64_t ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    double *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    double *RESTRICT _data_pdfs_20_34 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_312 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_33 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_313 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_35 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    double *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_32 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_317 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_31 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_318 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3;
    double *RESTRICT _data_force_20_31 =
        _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    double *RESTRICT _data_pdfs_20_316 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_311 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_36 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_314 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3;
    double *RESTRICT _data_force_20_32 =
        _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    double *RESTRICT _data_pdfs_20_315 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_20_37_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_37;
      double *RESTRICT _data_pdfs_20_310_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_310;
      double *RESTRICT _data_force_20_30_10 =
          _stride_force_1 * ctr_1 + _data_force_20_30;
      double *RESTRICT _data_pdfs_20_34_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      double *RESTRICT _data_pdfs_20_312_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_312;
      double *RESTRICT _data_pdfs_20_33_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      double *RESTRICT _data_pdfs_20_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_313;
      double *RESTRICT _data_pdfs_20_35_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_35;
      double *RESTRICT _data_pdfs_20_39_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_39;
      double *RESTRICT _data_pdfs_20_30_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      double *RESTRICT _data_pdfs_20_38_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_38;
      double *RESTRICT _data_pdfs_20_32_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_32;
      double *RESTRICT _data_pdfs_20_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_317;
      double *RESTRICT _data_pdfs_20_31_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_31;
      double *RESTRICT _data_pdfs_20_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_318;
      double *RESTRICT _data_force_20_31_10 =
          _stride_force_1 * ctr_1 + _data_force_20_31;
      double *RESTRICT _data_pdfs_20_316_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_316;
      double *RESTRICT _data_pdfs_20_311_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_311;
      double *RESTRICT _data_pdfs_20_36_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_36;
      double *RESTRICT _data_pdfs_20_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_314;
      double *RESTRICT _data_force_20_32_10 =
          _stride_force_1 * ctr_1 + _data_force_20_32;
      double *RESTRICT _data_pdfs_20_315_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_315;
      for (int64_t ctr_0 = 0; ctr_0 < _size_force_0; ctr_0 += 1) {
        const double xi_198 = _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0];
        const double xi_199 = _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0];
        const double xi_200 = _data_force_20_30_10[_stride_force_0 * ctr_0];
        const double xi_201 = _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0];
        const double xi_202 = _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0];
        const double xi_203 = _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0];
        const double xi_204 = _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0];
        const double xi_205 = _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0];
        const double xi_206 = _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0];
        const double xi_207 = _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0];
        const double xi_208 = _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0];
        const double xi_209 = _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0];
        const double xi_210 = _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0];
        const double xi_211 = _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0];
        const double xi_212 = _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0];
        const double xi_213 = _data_force_20_31_10[_stride_force_0 * ctr_0];
        const double xi_214 = _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0];
        const double xi_215 = _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0];
        const double xi_216 = _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0];
        const double xi_217 = _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0];
        const double xi_218 = _data_force_20_32_10[_stride_force_0 * ctr_0];
        const double xi_219 = _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0];
        const double xi_0 = xi_212 + xi_217;
        const double xi_1 = xi_0 + xi_201;
        const double xi_2 = xi_211 + xi_215 + xi_219;
        const double xi_3 = xi_202 + xi_205;
        const double xi_4 = xi_203 + xi_206;
        const double xi_5 = xi_209 + xi_214;
        const double xi_6 = xi_210 + xi_216;
        const double xi_8 = -xi_206;
        const double xi_9 = -xi_198 + xi_8;
        const double xi_10 = -xi_210;
        const double xi_11 = -xi_204;
        const double xi_12 = -xi_203;
        const double xi_13 = xi_10 + xi_11 + xi_12;
        const double xi_14 = -xi_209;
        const double xi_15 = -xi_199;
        const double xi_16 = xi_14 + xi_15;
        const double xi_17 = -xi_214;
        const double xi_18 = -xi_202;
        const double xi_19 = xi_17 + xi_18;
        const double xi_20 = -xi_212;
        const double xi_21 = xi_10 + xi_20;
        const double xi_22 = -xi_219;
        const double xi_23 = -xi_216;
        const double xi_24 = xi_17 + xi_215 + xi_22 + xi_23;
        const double xi_40 = xi_213 * 0.166666666666667;
        const double xi_48 = xi_200 * 0.166666666666667;
        const double xi_52 = xi_218 * 0.166666666666667;
        const double xi_55 = xi_213 * 0.5;
        const double xi_59 = xi_200 * 0.0833333333333333;
        const double xi_63 = xi_213 * 0.0833333333333333;
        const double xi_73 = xi_218 * 0.0833333333333333;
        const double xi_84 = -xi_207;
        const double xi_85 = xi_205 * 3.0 + xi_216 * 3.0 + xi_84;
        const double xi_86 =
            omega_even *
            (xi_202 * -3.0 + xi_209 * 3.0 + xi_211 * 3.0 + xi_214 * -3.0 +
             xi_215 * -3.0 + xi_219 * -3.0 + xi_85);
        const double xi_87 =
            xi_202 * 2.0 + xi_214 * 2.0 + xi_215 * 2.0 + xi_219 * 2.0;
        const double xi_88 = xi_201 * 5.0 + xi_203 * 5.0 + xi_87;
        const double xi_89 =
            omega_even *
            (xi_204 * -5.0 + xi_209 * -2.0 + xi_210 * -5.0 + xi_211 * -2.0 +
             xi_212 * -5.0 + xi_217 * -5.0 + xi_85 + xi_88);
        const double xi_92 = -xi_215;
        const double xi_93 = xi_18 + xi_92;
        const double xi_94 = -xi_208;
        const double xi_97 = -xi_217;
        const double xi_98 = xi_11 + xi_15 + xi_21 + xi_97;
        const double xi_100 = xi_204 * 2.0;
        const double xi_101 = xi_217 * 2.0;
        const double xi_102 = xi_210 * 2.0 + xi_212 * 2.0;
        const double xi_103 =
            omega_even *
            (xi_100 + xi_101 + xi_102 + xi_198 * -7.0 + xi_199 * -7.0 +
             xi_205 * -4.0 + xi_206 * -7.0 + xi_208 * -7.0 + xi_209 * 5.0 +
             xi_211 * 5.0 + xi_216 * -4.0 + xi_84 + xi_88);
        const double xi_104 = xi_202 + xi_92;
        const double xi_105 = xi_104 + xi_14 + xi_211 + xi_214 + xi_22;
        const double xi_107 = xi_105 * xi_106;
        const double xi_108 = xi_198 * 2.0;
        const double xi_109 = xi_199 * 2.0;
        const double xi_110 = xi_206 * 2.0 + xi_208 * -2.0;
        const double xi_111 = -xi_108 + xi_109 + xi_110 + xi_14 + xi_19 + xi_2;
        const double xi_113 = xi_111 * xi_112;
        const double xi_114 = -xi_113;
        const double xi_116 = xi_199 + xi_94;
        const double xi_120 = xi_204 + xi_210;
        const double xi_124 = xi_103 * -0.0198412698412698;
        const double xi_132 = xi_204 + xi_97;
        const double xi_133 = xi_12 + xi_132 + xi_20 + xi_201 + xi_210;
        const double xi_134 = xi_106 * xi_133;
        const double xi_135 = xi_1 + xi_108 - xi_109 + xi_110 + xi_13;
        const double xi_136 = xi_112 * xi_135;
        const double xi_138 = -xi_136;
        const double xi_139 = xi_214 + xi_219;
        const double xi_140 = xi_139 + xi_205 + xi_23 + xi_93;
        const double xi_141 = xi_106 * xi_140;
        const double xi_144 = -xi_100 - xi_101 + xi_102 + xi_24 + xi_3;
        const double xi_145 = xi_112 * xi_144;
        const double xi_146 = -xi_145;
        const double xi_148 = xi_145;
        const double xi_152 = xi_103 * 0.0138888888888889;
        const double xi_168 = xi_89 * -0.00714285714285714;
        const double xi_170 = xi_86 * 0.025;
        const double xi_173 = xi_144 * xi_172;
        const double xi_175 = xi_140 * xi_174;
        const double xi_176 = xi_103 * -0.00396825396825397;
        const double xi_180 = xi_111 * xi_172;
        const double xi_181 = xi_105 * xi_174;
        const double xi_187 = xi_89 * 0.0178571428571429;
        const double xi_190 = xi_133 * xi_174;
        const double xi_191 = xi_135 * xi_172;
        const double vel0Term = xi_1 + xi_199 + xi_208;
        const double vel1Term = xi_198 + xi_2;
        const double vel2Term = xi_204 + xi_3;
        const double rho =
            vel0Term + vel1Term + vel2Term + xi_207 + xi_4 + xi_5 + xi_6;
        const double xi_7 = 1 / (rho);
        const double u_0 = xi_7 * (vel0Term + xi_13 + xi_9);
        const double xi_25 = u_0 * xi_200;
        const double xi_26 = xi_25 * 0.333333333333333;
        const double xi_32 = -xi_26;
        const double xi_90 = rho * (u_0 * u_0);
        const double xi_129 = rho * u_0;
        const double xi_130 = -vel0Term + xi_120 + xi_129 + xi_198 + xi_4;
        const double xi_131 = xi_118 * xi_130;
        const double xi_158 = xi_130 * xi_154;
        const double u_1 = xi_7 * (vel1Term + xi_16 + xi_19 + xi_208 + xi_8);
        const double xi_27 = u_1 * xi_213;
        const double xi_28 = xi_27 * 0.333333333333333;
        const double xi_33 = -xi_28;
        const double xi_54 = u_1 * 0.5;
        const double xi_57 = xi_56 * (u_0 * xi_55 + xi_200 * xi_54);
        const double xi_58 = -xi_57;
        const double xi_95 = rho * (u_1 * u_1);
        const double xi_96 = xi_9 + xi_94 + xi_95;
        const double xi_115 = rho * u_1;
        const double xi_117 =
            -vel1Term + xi_115 + xi_116 + xi_202 + xi_206 + xi_5;
        const double xi_119 = xi_117 * xi_118;
        const double xi_150 = xi_149 * (u_0 * xi_115 + xi_116 + xi_198 + xi_8);
        const double xi_155 = xi_117 * xi_154;
        const double xi_156 = xi_155;
        const double xi_157 = xi_113 + xi_156;
        const double xi_166 = -xi_155;
        const double xi_167 = xi_114 + xi_166;
        const double xi_182 = xi_156 - xi_180 + xi_181;
        const double xi_183 = xi_166 + xi_180 - xi_181;
        const double u_2 = xi_7 * (vel2Term + xi_21 + xi_217 + xi_24);
        const double xi_29 = u_2 * xi_218;
        const double xi_30 = xi_29 * 0.333333333333333;
        const double xi_31 = (-omega_bulk + 2.0) * (xi_26 + xi_28 + xi_30);
        const double xi_34 = xi_29 * 0.666666666666667 + xi_32 + xi_33;
        const double xi_37 = -xi_30;
        const double xi_38 = xi_27 * 0.666666666666667 + xi_32 + xi_37;
        const double xi_39 = xi_25 * 0.666666666666667 + xi_33 + xi_37;
        const double xi_42 = xi_34 * xi_41;
        const double xi_43 = -xi_42;
        const double xi_44 = xi_39 * xi_41;
        const double xi_45 = -xi_44;
        const double xi_47 = xi_38 * xi_46 + xi_43 + xi_45;
        const double xi_49 = xi_38 * xi_41;
        const double xi_50 = -xi_49;
        const double xi_51 = xi_39 * xi_46 + xi_43 + xi_50;
        const double xi_53 = xi_34 * xi_46 + xi_45 + xi_50;
        const double xi_60 = xi_44 - xi_59;
        const double xi_62 = -xi_34 * xi_61;
        const double xi_64 = xi_31 * 0.125;
        const double xi_65 = xi_49 + xi_64;
        const double xi_66 = xi_63 + xi_65;
        const double xi_67 = xi_62 + xi_66;
        const double xi_68 = xi_44 + xi_59;
        const double xi_69 = -xi_63 + xi_65;
        const double xi_70 = xi_62 + xi_69;
        const double xi_71 = xi_56 * (u_2 * xi_55 + xi_218 * xi_54);
        const double xi_72 = -xi_39 * xi_61;
        const double xi_74 = xi_42 + xi_73;
        const double xi_75 = xi_72 + xi_74;
        const double xi_76 = -xi_71;
        const double xi_77 = xi_56 * (u_0 * xi_218 * 0.5 + u_2 * xi_200 * 0.5);
        const double xi_78 = -xi_77;
        const double xi_79 = -xi_38 * xi_61;
        const double xi_80 = xi_64 + xi_74 + xi_79;
        const double xi_81 = xi_42 - xi_73;
        const double xi_82 = xi_72 + xi_81;
        const double xi_83 = xi_64 + xi_79 + xi_81;
        const double xi_91 = rho * (u_2 * u_2);
        const double xi_99 = omega_bulk * (xi_17 + xi_207 + xi_22 + xi_90 +
                                           xi_91 + xi_93 + xi_96 + xi_98);
        const double xi_121 = xi_205 + xi_216 - xi_91;
        const double xi_122 =
            omega_shear * (xi_0 + xi_120 + xi_121 + xi_16 - xi_211 + xi_96);
        const double xi_123 = xi_122 * 0.125;
        const double xi_125 =
            omega_shear *
            (xi_121 + xi_201 * -2.0 + xi_203 * -2.0 + xi_209 + xi_211 + xi_87 +
             xi_9 + xi_90 * 2.0 + xi_94 - xi_95 + xi_98);
        const double xi_126 = xi_125 * -0.0416666666666667;
        const double xi_127 = xi_126 + xi_86 * -0.05;
        const double xi_128 =
            xi_123 + xi_124 + xi_127 + xi_89 * 0.0142857142857143;
        const double xi_137 =
            xi_124 + xi_125 * 0.0833333333333333 + xi_89 * -0.0357142857142857;
        const double xi_142 =
            rho * u_2 - vel2Term + xi_139 + xi_212 + xi_6 + xi_92 + xi_97;
        const double xi_143 = xi_118 * xi_142;
        const double xi_147 = xi_103 * 0.0158730158730159 - xi_123 + xi_127 +
                              xi_89 * -0.0214285714285714;
        const double xi_151 = xi_122 * 0.0625;
        const double xi_153 = -xi_150 + xi_151 + xi_152;
        const double xi_159 = xi_99 * 0.0416666666666667;
        const double xi_160 = xi_125 * 0.0208333333333333 + xi_159;
        const double xi_161 = -xi_158 + xi_160;
        const double xi_162 = xi_138 + xi_161;
        const double xi_163 = xi_150 + xi_151 + xi_152;
        const double xi_164 = xi_158 + xi_160;
        const double xi_165 = xi_136 + xi_164;
        const double xi_169 = xi_149 * (u_2 * xi_115 + xi_104 + xi_17 + xi_219);
        const double xi_171 = xi_126 + xi_159 + xi_168 + xi_169 + xi_170;
        const double xi_177 = xi_142 * xi_154;
        const double xi_178 = xi_176 + xi_177;
        const double xi_179 = -xi_173 + xi_175 + xi_178;
        const double xi_184 = xi_126 + xi_159 + xi_168 - xi_169 + xi_170;
        const double xi_185 = xi_149 * (u_2 * xi_129 + xi_10 + xi_132 + xi_212);
        const double xi_186 = -xi_151;
        const double xi_188 = -xi_185 + xi_186 + xi_187;
        const double xi_189 = xi_148 + xi_178;
        const double xi_192 = xi_161 - xi_190 + xi_191;
        const double xi_193 = xi_185 + xi_186 + xi_187;
        const double xi_194 = xi_164 + xi_190 - xi_191;
        const double xi_195 = xi_176 - xi_177;
        const double xi_196 = xi_173 - xi_175 + xi_195;
        const double xi_197 = xi_146 + xi_195;
        const double forceTerm_0 =
            xi_31 * -1.5 - xi_34 * xi_36 - xi_36 * xi_38 - xi_36 * xi_39;
        const double forceTerm_1 = xi_40 + xi_47;
        const double forceTerm_2 = -xi_40 + xi_47;
        const double forceTerm_3 = -xi_48 + xi_51;
        const double forceTerm_4 = xi_48 + xi_51;
        const double forceTerm_5 = xi_52 + xi_53;
        const double forceTerm_6 = -xi_52 + xi_53;
        const double forceTerm_7 = xi_58 + xi_60 + xi_67;
        const double forceTerm_8 = xi_57 + xi_67 + xi_68;
        const double forceTerm_9 = xi_57 + xi_60 + xi_70;
        const double forceTerm_10 = xi_58 + xi_68 + xi_70;
        const double forceTerm_11 = xi_66 + xi_71 + xi_75;
        const double forceTerm_12 = xi_69 + xi_75 + xi_76;
        const double forceTerm_13 = xi_60 + xi_78 + xi_80;
        const double forceTerm_14 = xi_68 + xi_77 + xi_80;
        const double forceTerm_15 = xi_66 + xi_76 + xi_82;
        const double forceTerm_16 = xi_69 + xi_71 + xi_82;
        const double forceTerm_17 = xi_60 + xi_77 + xi_83;
        const double forceTerm_18 = xi_68 + xi_78 + xi_83;
        _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_0 + xi_103 * 0.0238095238095238 + xi_207 + xi_86 * 0.1 +
            xi_89 * 0.0428571428571429 + xi_99 * -0.5;
        _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_1 - xi_107 + xi_114 + xi_119 + xi_128 + xi_211;
        _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_2 + xi_107 + xi_113 - xi_119 + xi_128 + xi_209;
        _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_3 - xi_131 + xi_134 + xi_136 + xi_137 + xi_203;
        _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_4 + xi_131 - xi_134 + xi_137 + xi_138 + xi_201;
        _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_5 - xi_141 + xi_143 + xi_146 + xi_147 + xi_205;
        _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_6 + xi_141 - xi_143 + xi_147 + xi_148 + xi_216;
        _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_7 + xi_153 + xi_157 + xi_162 + xi_198;
        _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_8 + xi_157 + xi_163 + xi_165 + xi_208;
        _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_9 + xi_162 + xi_163 + xi_167 + xi_206;
        _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_10 + xi_153 + xi_165 + xi_167 + xi_199;
        _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_11 + xi_171 + xi_179 + xi_182 + xi_215;
        _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_12 + xi_179 + xi_183 + xi_184 + xi_202;
        _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_13 + xi_188 + xi_189 + xi_192 + xi_204;
        _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_14 + xi_189 + xi_193 + xi_194 + xi_217;
        _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_15 + xi_182 + xi_184 + xi_196 + xi_219;
        _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_16 + xi_171 + xi_183 + xi_196 + xi_214;
        _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_17 + xi_192 + xi_193 + xi_197 + xi_210;
        _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0] =
            forceTerm_18 + xi_188 + xi_194 + xi_197 + xi_212;
      }
    }
  }
}
} // namespace internal_collidesweep

void CollideSweep::operator()(IBlock *block) {
  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);
  auto force = block->getData<field::GhostLayerField<double, 3>>(forceID);

  auto &omega_even = this->omega_even_;
  auto &omega_odd = this->omega_odd_;
  auto &omega_bulk = this->omega_bulk_;
  auto &omega_shear = this->omega_shear_;
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force->nrOfGhostLayers()));
  double *RESTRICT const _data_force = force->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()));
  double *RESTRICT _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
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
  internal_collidesweep::collidesweep(
      _data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2,
      _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3,
      _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3,
      omega_bulk, omega_even, omega_odd, omega_shear);
}

void CollideSweep::runOnCellInterval(
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

  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);
  auto force = block->getData<field::GhostLayerField<double, 3>>(forceID);

  auto &omega_even = this->omega_even_;
  auto &omega_odd = this->omega_odd_;
  auto &omega_bulk = this->omega_bulk_;
  auto &omega_shear = this->omega_shear_;
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(force->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(force->nrOfGhostLayers()));
  double *RESTRICT const _data_force =
      force->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(pdfs->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(pdfs->nrOfGhostLayers()));
  double *RESTRICT _data_pdfs =
      pdfs->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
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
  internal_collidesweep::collidesweep(
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