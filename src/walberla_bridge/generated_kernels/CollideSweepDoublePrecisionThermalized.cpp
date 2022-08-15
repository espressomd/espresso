// kernel generated with pystencils v1.0+12.g54b91e2, lbmpy
// v1.0+9.g19115d4.dirty, lbmpy_walberla/pystencils_walberla from commit
// e1fe2ad1dcbe8f31ea79d95e8a5a5cc0ee3691f3

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

#include <cmath>

#include "CollideSweepDoublePrecisionThermalized.h"
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

namespace internal_0d943397135d13b4628c5752888935d7 {
static FUNC_PREFIX void
collidesweepdoubleprecisionthermalized_collidesweepdoubleprecisionthermalized(
    double *RESTRICT const _data_force, double *RESTRICT _data_pdfs,
    int64_t const _size_force_0, int64_t const _size_force_1,
    int64_t const _size_force_2, int64_t const _stride_force_0,
    int64_t const _stride_force_1, int64_t const _stride_force_2,
    int64_t const _stride_force_3, int64_t const _stride_pdfs_0,
    int64_t const _stride_pdfs_1, int64_t const _stride_pdfs_2,
    int64_t const _stride_pdfs_3, uint32_t block_offset_0,
    uint32_t block_offset_1, uint32_t block_offset_2, double kT,
    double omega_bulk, double omega_even, double omega_odd, double omega_shear,
    double rho, uint32_t seed, uint32_t time_step) {
  const double xi_22 = omega_shear * -0.5 + 1.0;
  const double xi_59 = kT * rho;
  const double xi_60 = pow(
      xi_59 * (-1.0 * ((omega_even * -1.0 + 1.0) * (omega_even * -1.0 + 1.0)) +
               1.0),
      0.5);
  const double xi_63 = 2.4494897427831779;
  const bool xi_66 = kT > 0.0;
  const double xi_68 = pow(
      xi_59 *
          (-1.0 * ((omega_odd * -1.0 + 1.0) * (omega_odd * -1.0 + 1.0)) + 1.0),
      0.5);
  const double xi_69 = xi_68 * 1.4142135623730951;
  const double xi_70 = xi_69 * 0.5;
  const double xi_73 = xi_63 * xi_68;
  const double xi_74 = xi_73 * 0.16666666666666666;
  const double xi_79 = pow(xi_59 * (-1.0 * ((omega_shear * -1.0 + 1.0) *
                                            (omega_shear * -1.0 + 1.0)) +
                                    1.0),
                           0.5);
  const double xi_80 = xi_79 * 0.5;
  const double xi_111 = xi_73 * 0.083333333333333329;
  const double xi_113 = xi_69 * 0.25;
  for (int64_t ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1) {
    double *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_314 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_313 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_318 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_36 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    double *RESTRICT _data_force_20_31 =
        _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    double *RESTRICT _data_pdfs_20_32 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_317 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_315 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_31 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    double *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    double *RESTRICT _data_pdfs_20_312 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_35 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_33 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_316 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    double *RESTRICT _data_force_20_32 =
        _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    double *RESTRICT _data_pdfs_20_311 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3;
    double *RESTRICT _data_pdfs_20_34 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1) {
      double *RESTRICT _data_pdfs_20_37_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_37;
      double *RESTRICT _data_pdfs_20_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_314;
      double *RESTRICT _data_pdfs_20_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_313;
      double *RESTRICT _data_pdfs_20_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_318;
      double *RESTRICT _data_pdfs_20_36_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_36;
      double *RESTRICT _data_pdfs_20_30_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      double *RESTRICT _data_force_20_31_10 =
          _stride_force_1 * ctr_1 + _data_force_20_31;
      double *RESTRICT _data_pdfs_20_32_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_32;
      double *RESTRICT _data_pdfs_20_38_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_38;
      double *RESTRICT _data_pdfs_20_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_317;
      double *RESTRICT _data_pdfs_20_315_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_315;
      double *RESTRICT _data_pdfs_20_39_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_39;
      double *RESTRICT _data_pdfs_20_31_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_31;
      double *RESTRICT _data_force_20_30_10 =
          _stride_force_1 * ctr_1 + _data_force_20_30;
      double *RESTRICT _data_pdfs_20_312_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_312;
      double *RESTRICT _data_pdfs_20_35_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_35;
      double *RESTRICT _data_pdfs_20_33_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      double *RESTRICT _data_pdfs_20_316_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_316;
      double *RESTRICT _data_pdfs_20_310_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_310;
      double *RESTRICT _data_force_20_32_10 =
          _stride_force_1 * ctr_1 + _data_force_20_32;
      double *RESTRICT _data_pdfs_20_311_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_311;
      double *RESTRICT _data_pdfs_20_34_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      for (int64_t ctr_0 = 0; ctr_0 < _size_force_0; ctr_0 += 1) {
        const double xi_132 = _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0];
        const double xi_133 = _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0];
        const double xi_134 = _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0];
        const double xi_135 = _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0];
        const double xi_136 = _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0];
        const double xi_137 = _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0];
        const double xi_138 = _data_force_20_31_10[_stride_force_0 * ctr_0];
        const double xi_139 = _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0];
        const double xi_140 = _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0];
        const double xi_141 = _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0];
        const double xi_142 = _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0];
        const double xi_143 = _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0];
        const double xi_144 = _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0];
        const double xi_145 = _data_force_20_30_10[_stride_force_0 * ctr_0];
        const double xi_146 = _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0];
        const double xi_147 = _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0];
        const double xi_148 = _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0];
        const double xi_149 = _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0];
        const double xi_150 = _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0];
        const double xi_151 = _data_force_20_32_10[_stride_force_0 * ctr_0];
        const double xi_152 = _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0];
        const double xi_153 = _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0];

        double random_7_0{};
        double random_7_1{};
        if (kT > 0.) {
          philox_double2(time_step, block_offset_0 + ctr_0,
                         block_offset_1 + ctr_1, block_offset_2 + ctr_2, 7,
                         seed, random_7_0, random_7_1);
        }

        double random_6_0{};
        double random_6_1{};
        if (kT > 0.) {
          philox_double2(time_step, block_offset_0 + ctr_0,
                         block_offset_1 + ctr_1, block_offset_2 + ctr_2, 6,
                         seed, random_6_0, random_6_1);
        }

        double random_5_0{};
        double random_5_1{};
        if (kT > 0.) {
          philox_double2(time_step, block_offset_0 + ctr_0,
                         block_offset_1 + ctr_1, block_offset_2 + ctr_2, 5,
                         seed, random_5_0, random_5_1);
        }

        double random_4_0{};
        double random_4_1{};
        if (kT > 0.) {
          philox_double2(time_step, block_offset_0 + ctr_0,
                         block_offset_1 + ctr_1, block_offset_2 + ctr_2, 4,
                         seed, random_4_0, random_4_1);
        }

        double random_3_0{};
        double random_3_1{};
        if (kT > 0.) {
          philox_double2(time_step, block_offset_0 + ctr_0,
                         block_offset_1 + ctr_1, block_offset_2 + ctr_2, 3,
                         seed, random_3_0, random_3_1);
        }

        double random_2_0{};
        double random_2_1{};
        if (kT > 0.) {
          philox_double2(time_step, block_offset_0 + ctr_0,
                         block_offset_1 + ctr_1, block_offset_2 + ctr_2, 2,
                         seed, random_2_0, random_2_1);
        }

        double random_1_0{};
        double random_1_1{};
        if (kT > 0.) {
          philox_double2(time_step, block_offset_0 + ctr_0,
                         block_offset_1 + ctr_1, block_offset_2 + ctr_2, 1,
                         seed, random_1_0, random_1_1);
        }

        double random_0_0{};
        double random_0_1{};
        if (kT > 0.) {
          philox_double2(time_step, block_offset_0 + ctr_0,
                         block_offset_1 + ctr_1, block_offset_2 + ctr_2, 0,
                         seed, random_0_0, random_0_1);
        }
        const double xi_0 = xi_134 + xi_141;
        const double xi_1 = xi_132 + xi_143;
        const double xi_2 = xi_146 + xi_149;
        const double xi_3 = xi_136 + xi_147;
        const double xi_4 = xi_142 + xi_152;
        const double xi_6 = xi_133 + xi_135;
        const double xi_7 = xi_140 + xi_150;
        const double xi_61 = xi_60 * (random_6_0 - 0.5) * 3.7416573867739413;
        const double xi_62 = xi_60 * (random_7_0 - 0.5) * 5.4772255750516612;
        const double xi_64 = xi_63 * (random_2_1 - 0.5) *
                             pow(xi_59 * (-1.0 * ((omega_bulk * -1.0 + 1.0) *
                                                  (omega_bulk * -1.0 + 1.0)) +
                                          1.0),
                                 0.5);
        const double xi_65 = xi_60 * (random_6_1 - 0.5) * 8.3666002653407556;
        const double xi_67 = random_5_1 - 0.5;
        const double xi_71 = xi_67 * xi_70;
        const double xi_72 = random_3_0 - 0.5;
        const double xi_75 = xi_72 * xi_74;
        const double xi_76 = xi_75 * -1.0;
        const double xi_77 = xi_61 * -0.11904761904761904;
        const double xi_78 = random_0_1 - 0.5;
        const double xi_81 = xi_78 * xi_80;
        const double xi_82 = xi_79 * (random_0_0 - 0.5) * 1.7320508075688772;
        const double xi_83 = xi_82 * -0.16666666666666666;
        const double xi_84 = xi_62 * -0.10000000000000001 + xi_83;
        const double xi_85 =
            xi_65 * 0.028571428571428571 + xi_77 + xi_81 + xi_84;
        const double xi_86 = random_4_1 - 0.5;
        const double xi_87 = xi_70 * xi_86;
        const double xi_88 = random_4_0 - 0.5;
        const double xi_89 = xi_74 * xi_88;
        const double xi_90 =
            xi_65 * -0.071428571428571425 + xi_77 + xi_82 * 0.33333333333333331;
        const double xi_91 = xi_89 * -1.0;
        const double xi_92 = random_5_0 - 0.5;
        const double xi_93 = xi_70 * xi_92;
        const double xi_94 = random_3_1 - 0.5;
        const double xi_95 = xi_74 * xi_94;
        const double xi_96 = xi_95 * -1.0;
        const double xi_97 = xi_61 * 0.095238095238095233 +
                             xi_65 * -0.042857142857142858 + xi_81 * -1.0 +
                             xi_84;
        const double xi_98 = xi_80 * (random_1_0 - 0.5);
        const double xi_99 = xi_98 * -1.0;
        const double xi_100 = xi_78 * xi_79 * 0.25;
        const double xi_101 = xi_61 * 0.083333333333333329;
        const double xi_102 = xi_64 * 0.083333333333333329;
        const double xi_103 = xi_102 + xi_82 * 0.083333333333333329;
        const double xi_104 = xi_100 + xi_101 + xi_103 + xi_75;
        const double xi_105 = xi_100 + xi_101 + xi_103 + xi_76;
        const double xi_106 = xi_61 * -0.023809523809523808;
        const double xi_107 = xi_80 * (random_2_0 - 0.5);
        const double xi_108 = xi_65 * -0.014285714285714285;
        const double xi_109 = xi_62 * 0.050000000000000003;
        const double xi_110 =
            xi_102 + xi_106 + xi_107 + xi_108 + xi_109 + xi_83;
        const double xi_112 = xi_111 * xi_94;
        const double xi_114 = xi_113 * xi_92;
        const double xi_115 = xi_112 * -1.0 + xi_114;
        const double xi_116 = xi_111 * xi_72;
        const double xi_117 = xi_113 * xi_67;
        const double xi_118 = xi_116 * -1.0 + xi_117;
        const double xi_119 = xi_116 + xi_117 * -1.0;
        const double xi_120 =
            xi_102 + xi_106 + xi_107 * -1.0 + xi_108 + xi_109 + xi_83;
        const double xi_121 = xi_80 * (random_1_1 - 0.5);
        const double xi_122 = xi_121 * -1.0;
        const double xi_123 = xi_100 * -1.0;
        const double xi_124 = xi_65 * 0.035714285714285712;
        const double xi_125 = xi_103 + xi_106 + xi_123 + xi_124 + xi_95;
        const double xi_126 = xi_113 * xi_86;
        const double xi_127 = xi_111 * xi_88;
        const double xi_128 = xi_126 * -1.0 + xi_127;
        const double xi_129 = xi_126 + xi_127 * -1.0;
        const double xi_130 = xi_112 + xi_114 * -1.0;
        const double xi_131 = xi_103 + xi_106 + xi_123 + xi_124 + xi_96;
        const double partial_m_m1_0_e_0 = xi_0 + xi_148;
        const double partial_m_m1_e_0_0 = partial_m_m1_0_e_0 + xi_1;
        const double partial_m_0_m1_e_0 = xi_139 + xi_2;
        const double partial_m_0_0_e_0 = xi_137 + xi_3;
        const double partial_m_0_1_e_0 = xi_144 + xi_4;
        const double xi_5 = partial_m_0_1_e_0 + partial_m_0_m1_e_0;
        const double partial_m_0_e_0_0 = partial_m_0_0_e_0 + xi_5;
        const double partial_m_1_0_e_0 = xi_153 + xi_6;
        const double partial_m_1_e_0_0 = partial_m_1_0_e_0 + xi_7;
        const double xi_10 = partial_m_1_e_0_0 + partial_m_m1_e_0_0;
        const double partial_m_m1_e_1_0 = xi_132 + xi_143 * -1.0;
        const double partial_m_0_e_1_0 =
            partial_m_0_1_e_0 + partial_m_0_m1_e_0 * -1.0;
        const double partial_m_1_e_1_0 = xi_140 + xi_150 * -1.0;
        const double xi_11 = partial_m_1_e_1_0 + partial_m_m1_e_1_0;
        const double partial_m_m1_0_e_1 = xi_134 + xi_141 * -1.0;
        const double partial_m_0_m1_e_1 = xi_146 + xi_149 * -1.0;
        const double partial_m_0_0_e_1 = xi_136 * -1.0 + xi_147;
        const double partial_m_0_1_e_1 = xi_142 * -1.0 + xi_152;
        const double xi_8 = partial_m_0_1_e_1 + partial_m_0_m1_e_1;
        const double partial_m_0_e_0_1 = partial_m_0_0_e_1 + xi_8;
        const double partial_m_1_0_e_1 = xi_133 + xi_135 * -1.0;
        const double xi_12 = partial_m_1_0_e_1 + partial_m_m1_0_e_1;
        const double partial_m_m1_e_2_0 = xi_1;
        const double partial_m_0_e_2_0 = xi_5;
        const double partial_m_1_e_2_0 = xi_7;
        const double xi_13 = partial_m_1_e_2_0 + partial_m_m1_e_2_0;
        const double partial_m_m1_0_e_2 = xi_0;
        const double partial_m_0_m1_e_2 = xi_2;
        const double partial_m_0_0_e_2 = xi_3;
        const double partial_m_0_1_e_2 = xi_4;
        const double xi_9 = partial_m_0_1_e_2 + partial_m_0_m1_e_2;
        const double partial_m_0_e_0_2 = partial_m_0_0_e_2 + xi_9;
        const double partial_m_1_0_e_2 = xi_6;
        const double xi_14 = partial_m_1_0_e_2 + partial_m_m1_0_e_2;
        const double partial_m_0_e_1_1 =
            partial_m_0_1_e_1 + partial_m_0_m1_e_1 * -1.0;
        const double partial_m_0_e_2_1 = xi_8;
        const double partial_m_0_e_1_2 =
            partial_m_0_1_e_2 + partial_m_0_m1_e_2 * -1.0;
        const double partial_m_0_e_2_2 = xi_9;
        const double m_000 = partial_m_0_e_0_0 + xi_10;
        const double m_100 = partial_m_1_e_0_0 + partial_m_m1_e_0_0 * -1.0;
        const double xi_19 = m_100 * -1.0;
        const double m_010 = partial_m_0_e_1_0 + xi_11;
        const double xi_17 = m_010 * -1.0;
        const double m_001 = partial_m_0_e_0_1 + xi_12;
        const double xi_18 = m_001 * -1.0;
        const double m_200 = xi_10;
        const double xi_20 = m_000 + m_200 * -6.0;
        const double m_020 = partial_m_0_e_2_0 + xi_13;
        const double m_002 = partial_m_0_e_0_2 + xi_14;
        const double xi_16 = m_002 * -1.0;
        const double m_110 = partial_m_1_e_1_0 + partial_m_m1_e_1_0 * -1.0;
        const double m_101 = partial_m_1_0_e_1 + partial_m_m1_0_e_1 * -1.0;
        const double m_210 = xi_11;
        const double m_201 = xi_12;
        const double m_120 = partial_m_1_e_2_0 + partial_m_m1_e_2_0 * -1.0;
        const double m_102 = partial_m_1_0_e_2 + partial_m_m1_0_e_2 * -1.0;
        const double m_220 = xi_13;
        const double xi_21 = m_002 * -4.0 + m_220 * 4.0;
        const double m_202 = xi_14;
        const double sub_f_to_m_0 = ((1.0) / (m_000));
        const double xi_15 = sub_f_to_m_0 * 0.5;
        const double u_0 = m_100 * sub_f_to_m_0 + xi_145 * xi_15;
        const double xi_23 = u_0 * xi_145;
        const double xi_27 = m_000 * (u_0 * u_0);
        const double xi_31 = m_000 * u_0;
        const double u_1 = m_010 * sub_f_to_m_0 + xi_138 * xi_15;
        const double xi_24 = u_1 * xi_138 * 2.0;
        const double xi_28 = m_000 * (u_1 * u_1);
        const double u_2 = m_001 * sub_f_to_m_0 + xi_15 * xi_151;
        const double xi_25 = u_2 * xi_151 * 2.0;
        const double xi_26 = xi_25 * -1.0;
        const double xi_29 = m_000 * (u_2 * u_2);
        const double xi_30 = xi_29 * -1.0;
        const double xi_32 = xi_29 * 0.66666666666666663;
        const double M_4 = m_020 * -1.0 + m_200 * 2.0 + xi_16;
        const double M_5 = m_020 + xi_16;
        const double M_9 = m_000 * -1.0 + m_002 + m_020 + m_200;
        const double M_10 = m_210 * 3.0 + xi_17;
        const double M_11 = m_201 * 3.0 + xi_18;
        const double M_12 = m_120 * 3.0 + xi_19;
        const double M_13 = m_102 * 2.0 + m_120 + xi_19;
        const double M_14 = m_201 + partial_m_0_e_2_1 * 2.0 + xi_18;
        const double M_15 = m_210 + partial_m_0_e_1_2 * 2.0 + xi_17;
        const double M_16 = m_002 * 3.0 + m_020 * -6.0 + m_220 * 18.0 + xi_20;
        const double M_17 = m_020 + m_202 * 14.0 + xi_20 + xi_21;
        const double M_18 = m_000 + m_020 * -4.0 + m_200 * -1.0 + m_202 * 4.0 +
                            partial_m_0_e_2_2 * 10.0 + xi_21;
        const double M_post_1 = m_100 + xi_145;
        const double xi_39 = M_post_1 * 0.33333333333333331;
        const double M_post_2 = m_010 + xi_138;
        const double xi_37 = M_post_2 * 0.33333333333333331;
        const double M_post_3 = m_001 + xi_151;
        const double xi_38 = M_post_3 * 0.33333333333333331;
        const double M_post_4 =
            M_4 +
            omega_shear * (M_4 * -1.0 + xi_27 * 2.0 + xi_28 * -1.0 + xi_30) +
            xi_22 * (xi_23 * 4.0 + xi_24 * -1.0 + xi_26);
        const double xi_35 = M_post_4 * -0.16666666666666666;
        const double M_post_5 = M_5 +
                                omega_shear * (M_5 * -1.0 + xi_28 + xi_30) +
                                xi_22 * (xi_24 + xi_26);
        const double xi_34 = M_post_5 * 0.5;
        const double xi_40 = M_post_5 * 0.25;
        const double M_post_6 = m_110 +
                                omega_shear * (m_110 * -1.0 + u_1 * xi_31) +
                                xi_22 * (u_0 * xi_138 + u_1 * xi_145);
        const double xi_47 = M_post_6 * 0.25;
        const double M_post_7 = m_101 +
                                omega_shear * (m_101 * -1.0 + u_2 * xi_31) +
                                xi_22 * (u_0 * xi_151 + u_2 * xi_145);
        const double xi_55 = M_post_7 * 0.25;
        const double M_post_8 =
            omega_shear * (m_000 * u_1 * u_2 + partial_m_0_e_1_1 * -1.0) +
            partial_m_0_e_1_1 + xi_22 * (u_1 * xi_151 + u_2 * xi_138);
        const double xi_51 = M_post_8 * 0.25;
        const double M_post_9 =
            M_9 + omega_bulk * (M_9 * -1.0 + xi_27 + xi_28 + xi_29) +
            (omega_bulk * -0.5 + 1.0) * (xi_23 * 2.0 + xi_24 + xi_25);
        const double xi_33 =
            M_post_9 * 0.33333333333333331 + m_000 * 0.33333333333333331;
        const double xi_36 = xi_33 + xi_35;
        const double xi_41 =
            M_post_9 * 0.16666666666666666 + m_000 * 0.1111111111111111;
        const double xi_42 = M_post_4 * 0.083333333333333329 + xi_41;
        const double M_post_10 = M_10 * omega_odd * -1.0 + M_10;
        const double M_post_11 = M_11 * omega_odd * -1.0 + M_11;
        const double M_post_12 = M_12 * omega_odd * -1.0 + M_12;
        const double M_post_13 = M_13 * omega_odd * -1.0 + M_13;
        const double M_post_14 = M_14 * omega_odd * -1.0 + M_14;
        const double M_post_15 = M_15 * omega_odd * -1.0 + M_15;
        const double M_post_16 =
            M_16 + omega_even * (M_16 * -1.0 + xi_29 * 3.0);
        const double xi_43 = M_post_16 * -0.015873015873015872;
        const double M_post_17 =
            M_17 +
            omega_even * (M_17 * -1.0 + xi_28 * 2.3333333333333335 + xi_32);
        const double M_post_18 =
            M_18 + omega_even * (M_18 * -1.0 + xi_27 * 1.6666666666666667 +
                                 xi_28 * 0.66666666666666663 + xi_32);
        const double m_post_200 = M_post_4 * 0.33333333333333331 + xi_33;
        const double m_post_020 = xi_34 + xi_36;
        const double m_post_002 = xi_34 * -1.0 + xi_36;
        const double m_post_210 = M_post_10 * 0.33333333333333331 + xi_37;
        const double xi_50 = m_post_210 * 0.25;
        const double m_post_201 = M_post_11 * 0.33333333333333331 + xi_38;
        const double xi_58 = m_post_201 * 0.25;
        const double m_post_120 = M_post_12 * 0.33333333333333331 + xi_39;
        const double xi_49 = m_post_120 * 0.25;
        const double m_post_102 =
            M_post_12 * -0.16666666666666666 + M_post_13 * 0.5 + xi_39;
        const double xi_57 = m_post_102 * 0.25;
        const double m_post_021 =
            M_post_11 * -0.16666666666666666 + M_post_14 * 0.5 + xi_38;
        const double xi_54 = m_post_021 * 0.25;
        const double m_post_012 =
            M_post_10 * -0.16666666666666666 + M_post_15 * 0.5 + xi_37;
        const double xi_53 = m_post_012 * 0.25;
        const double m_post_220 =
            M_post_16 * 0.055555555555555552 + xi_40 + xi_42;
        const double xi_45 = m_post_220 * -0.5;
        const double xi_48 = m_post_220 * 0.25;
        const double m_post_202 =
            M_post_17 * 0.071428571428571425 + xi_40 * -1.0 + xi_42 + xi_43;
        const double xi_46 = m_post_202 * -0.5;
        const double xi_56 = m_post_202 * 0.25;
        const double m_post_022 = M_post_17 * -0.028571428571428571 +
                                  M_post_18 * 0.10000000000000001 + xi_35 +
                                  xi_41 + xi_43;
        const double xi_44 = m_post_022 * -0.5;
        const double xi_52 = m_post_022 * 0.25;
        const double sub_k_to_f_20 = m_post_020 * 0.5 + xi_44 + xi_45;
        const double sub_k_to_f_21 =
            M_post_2 * 0.5 + m_post_012 * -0.5 + m_post_210 * -0.5;
        const double sub_k_to_f_22 = m_post_200 * 0.5 + xi_45 + xi_46;
        const double sub_k_to_f_23 =
            M_post_1 * -0.5 + m_post_102 * 0.5 + m_post_120 * 0.5;
        const double sub_k_to_f_24 = m_post_002 * 0.5 + xi_44 + xi_46;
        const double sub_k_to_f_25 =
            M_post_3 * 0.5 + m_post_021 * -0.5 + m_post_201 * -0.5;
        const double sub_k_to_f_26 = xi_47 * -1.0 + xi_48;
        const double sub_k_to_f_27 = xi_49 * -1.0 + xi_50;
        const double sub_k_to_f_28 = xi_47 + xi_48;
        const double sub_k_to_f_29 = xi_49 + xi_50;
        const double sub_k_to_f_30 = xi_51 + xi_52;
        const double sub_k_to_f_31 = xi_53 + xi_54;
        const double sub_k_to_f_32 = xi_51 * -1.0 + xi_52;
        const double sub_k_to_f_33 = xi_53 * -1.0 + xi_54;
        const double sub_k_to_f_34 = xi_55 * -1.0 + xi_56;
        const double sub_k_to_f_35 = xi_57 * -1.0 + xi_58;
        const double sub_k_to_f_36 = xi_55 + xi_56;
        const double sub_k_to_f_37 = xi_57 + xi_58;
        _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0] =
            m_000 + m_post_002 * -1.0 + m_post_020 * -1.0 + m_post_022 +
            m_post_200 * -1.0 + m_post_202 + m_post_220 +
            ((xi_66)
                 ? (xi_61 * 0.14285714285714285 + xi_62 * 0.20000000000000001 +
                    xi_64 * -1.0 + xi_65 * 0.085714285714285715)
                 : (0.0));
        _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_20 + sub_k_to_f_21 +
            ((xi_66) ? (xi_71 * -1.0 + xi_76 + xi_85) : (0.0));
        _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_20 + sub_k_to_f_21 * -1.0 +
            ((xi_66) ? (xi_71 + xi_75 + xi_85) : (0.0));
        _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_22 + sub_k_to_f_23 +
            ((xi_66) ? (xi_87 + xi_89 + xi_90) : (0.0));
        _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_22 + sub_k_to_f_23 * -1.0 +
            ((xi_66) ? (xi_87 * -1.0 + xi_90 + xi_91) : (0.0));
        _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_24 + sub_k_to_f_25 +
            ((xi_66) ? (xi_93 * -1.0 + xi_96 + xi_97) : (0.0));
        _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_24 + sub_k_to_f_25 * -1.0 +
            ((xi_66) ? (xi_93 + xi_95 + xi_97) : (0.0));
        _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_26 + sub_k_to_f_27 +
            ((xi_66) ? (xi_104 + xi_91 + xi_99) : (0.0));
        _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_28 + sub_k_to_f_29 +
            ((xi_66) ? (xi_104 + xi_89 + xi_98) : (0.0));
        _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_28 + sub_k_to_f_29 * -1.0 +
            ((xi_66) ? (xi_105 + xi_91 + xi_98) : (0.0));
        _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_26 + sub_k_to_f_27 * -1.0 +
            ((xi_66) ? (xi_105 + xi_89 + xi_99) : (0.0));
        _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_30 + sub_k_to_f_31 +
            ((xi_66) ? (xi_110 + xi_115 + xi_118) : (0.0));
        _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_32 + sub_k_to_f_33 +
            ((xi_66) ? (xi_115 + xi_119 + xi_120) : (0.0));
        _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_34 + sub_k_to_f_35 +
            ((xi_66) ? (xi_122 + xi_125 + xi_128) : (0.0));
        _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_36 + sub_k_to_f_37 +
            ((xi_66) ? (xi_121 + xi_125 + xi_129) : (0.0));
        _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_32 + sub_k_to_f_33 * -1.0 +
            ((xi_66) ? (xi_118 + xi_120 + xi_130) : (0.0));
        _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_30 + sub_k_to_f_31 * -1.0 +
            ((xi_66) ? (xi_110 + xi_119 + xi_130) : (0.0));
        _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_36 + sub_k_to_f_37 * -1.0 +
            ((xi_66) ? (xi_121 + xi_128 + xi_131) : (0.0));
        _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0] =
            sub_k_to_f_34 + sub_k_to_f_35 * -1.0 +
            ((xi_66) ? (xi_122 + xi_129 + xi_131) : (0.0));
      }
    }
  }
}
} // namespace internal_0d943397135d13b4628c5752888935d7

void CollideSweepDoublePrecisionThermalized::run(IBlock *block) {
  auto force = block->getData<field::GhostLayerField<double, 3>>(forceID);
  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);

  auto &rho = this->rho_;
  auto &omega_odd = this->omega_odd_;
  auto block_offset_0 = this->block_offset_0_;
  auto &kT = this->kT_;
  auto block_offset_2 = this->block_offset_2_;
  auto &omega_bulk = this->omega_bulk_;
  auto &omega_shear = this->omega_shear_;
  auto &omega_even = this->omega_even_;
  auto &time_step = this->time_step_;
  auto block_offset_1 = this->block_offset_1_;
  auto &seed = this->seed_;
  block_offset_generator(block, block_offset_0, block_offset_1, block_offset_2);
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
  internal_0d943397135d13b4628c5752888935d7::
      collidesweepdoubleprecisionthermalized_collidesweepdoubleprecisionthermalized(
          _data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2,
          _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3,
          _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3,
          block_offset_0, block_offset_1, block_offset_2, kT, omega_bulk,
          omega_even, omega_odd, omega_shear, rho, seed, time_step);
}

void CollideSweepDoublePrecisionThermalized::runOnCellInterval(
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

  auto force = block->getData<field::GhostLayerField<double, 3>>(forceID);
  auto pdfs = block->getData<field::GhostLayerField<double, 19>>(pdfsID);

  auto &rho = this->rho_;
  auto &omega_odd = this->omega_odd_;
  auto block_offset_0 = this->block_offset_0_;
  auto &kT = this->kT_;
  auto block_offset_2 = this->block_offset_2_;
  auto &omega_bulk = this->omega_bulk_;
  auto &omega_shear = this->omega_shear_;
  auto &omega_even = this->omega_even_;
  auto &time_step = this->time_step_;
  auto block_offset_1 = this->block_offset_1_;
  auto &seed = this->seed_;
  block_offset_generator(block, block_offset_0, block_offset_1, block_offset_2);
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
  internal_0d943397135d13b4628c5752888935d7::
      collidesweepdoubleprecisionthermalized_collidesweepdoubleprecisionthermalized(
          _data_force, _data_pdfs, _size_force_0, _size_force_1, _size_force_2,
          _stride_force_0, _stride_force_1, _stride_force_2, _stride_force_3,
          _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2, _stride_pdfs_3,
          block_offset_0, block_offset_1, block_offset_2, kT, omega_bulk,
          omega_even, omega_odd, omega_shear, rho, seed, time_step);
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