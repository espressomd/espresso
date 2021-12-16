// kernel generated with pystencils v0.4.3+4.g30da657, lbmpy v0.4.3+2.g0e17e61,
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
//! \\file CollideSweepSinglePrecisionLeesEdwards.cpp
//! \\ingroup lbm
//! \\author lbmpy
//======================================================================================================================

#include <cmath>

#include "CollideSweepSinglePrecisionLeesEdwards.h"
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

namespace internal_ab1f3bc3368574afb482da84ccb58898 {
static FUNC_PREFIX void
collidesweepsingleprecisionleesedwards_collidesweepsingleprecisionleesedwards(
    float *RESTRICT const _data_force, float *RESTRICT _data_pdfs,
    float *RESTRICT const _data_velocity, int64_t const _size_force_0,
    int64_t const _size_force_1, int64_t const _size_force_2,
    int64_t const _stride_force_0, int64_t const _stride_force_1,
    int64_t const _stride_force_2, int64_t const _stride_force_3,
    int64_t const _stride_pdfs_0, int64_t const _stride_pdfs_1,
    int64_t const _stride_pdfs_2, int64_t const _stride_pdfs_3,
    int64_t const _stride_velocity_0, int64_t const _stride_velocity_1,
    int64_t const _stride_velocity_2, int64_t const _stride_velocity_3,
    float omega_shear) {
  for (int64_t ctr_2 = 0; ctr_2 < _size_force_2; ctr_2 += 1) {
    float *RESTRICT _data_pdfs_20_34 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 4 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_318 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 18 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_35 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 5 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_313 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 13 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_310 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 10 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_31 =
        _data_force + _stride_force_2 * ctr_2 + _stride_force_3;
    float *RESTRICT _data_velocity_20_32 =
        _data_velocity + _stride_velocity_2 * ctr_2 + 2 * _stride_velocity_3;
    float *RESTRICT _data_pdfs_20_316 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 16 * _stride_pdfs_3;
    float *RESTRICT _data_velocity_20_31 =
        _data_velocity + _stride_velocity_2 * ctr_2 + _stride_velocity_3;
    float *RESTRICT _data_velocity_20_30 =
        _data_velocity + _stride_velocity_2 * ctr_2;
    float *RESTRICT _data_pdfs_20_317 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 17 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_38 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 8 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_36 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 6 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_37 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 7 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_32 =
        _data_force + _stride_force_2 * ctr_2 + 2 * _stride_force_3;
    float *RESTRICT _data_pdfs_20_312 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 12 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_33 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 3 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_30 = _data_pdfs + _stride_pdfs_2 * ctr_2;
    float *RESTRICT _data_pdfs_20_311 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 11 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_315 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 15 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_39 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 9 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_314 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 14 * _stride_pdfs_3;
    float *RESTRICT _data_force_20_30 = _data_force + _stride_force_2 * ctr_2;
    float *RESTRICT _data_pdfs_20_32 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + 2 * _stride_pdfs_3;
    float *RESTRICT _data_pdfs_20_31 =
        _data_pdfs + _stride_pdfs_2 * ctr_2 + _stride_pdfs_3;
    for (int64_t ctr_1 = 0; ctr_1 < _size_force_1; ctr_1 += 1) {
      float *RESTRICT _data_pdfs_20_34_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_34;
      float *RESTRICT _data_pdfs_20_318_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_318;
      float *RESTRICT _data_pdfs_20_35_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_35;
      float *RESTRICT _data_pdfs_20_313_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_313;
      float *RESTRICT _data_pdfs_20_310_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_310;
      float *RESTRICT _data_force_20_31_10 =
          _stride_force_1 * ctr_1 + _data_force_20_31;
      float *RESTRICT _data_velocity_20_32_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_32;
      float *RESTRICT _data_pdfs_20_316_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_316;
      float *RESTRICT _data_velocity_20_31_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_31;
      float *RESTRICT _data_velocity_20_30_10 =
          _stride_velocity_1 * ctr_1 + _data_velocity_20_30;
      float *RESTRICT _data_pdfs_20_317_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_317;
      float *RESTRICT _data_pdfs_20_38_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_38;
      float *RESTRICT _data_pdfs_20_36_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_36;
      float *RESTRICT _data_pdfs_20_37_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_37;
      float *RESTRICT _data_force_20_32_10 =
          _stride_force_1 * ctr_1 + _data_force_20_32;
      float *RESTRICT _data_pdfs_20_312_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_312;
      float *RESTRICT _data_pdfs_20_33_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_33;
      float *RESTRICT _data_pdfs_20_30_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_30;
      float *RESTRICT _data_pdfs_20_311_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_311;
      float *RESTRICT _data_pdfs_20_315_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_315;
      float *RESTRICT _data_pdfs_20_39_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_39;
      float *RESTRICT _data_pdfs_20_314_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_314;
      float *RESTRICT _data_force_20_30_10 =
          _stride_force_1 * ctr_1 + _data_force_20_30;
      float *RESTRICT _data_pdfs_20_32_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_32;
      float *RESTRICT _data_pdfs_20_31_10 =
          _stride_pdfs_1 * ctr_1 + _data_pdfs_20_31;
      for (int64_t ctr_0 = 0; ctr_0 < _size_force_0; ctr_0 += 1) {
        const float xi_26 = _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0];
        const float xi_27 = _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0];
        const float xi_28 = _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0];
        const float xi_29 = _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0];
        const float xi_30 = _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0];
        const float xi_31 = _data_force_20_31_10[_stride_force_0 * ctr_0];
        const float xi_32 = _data_velocity_20_32_10[_stride_velocity_0 * ctr_0];
        const float xi_33 = _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0];
        const float xi_34 = _data_velocity_20_31_10[_stride_velocity_0 * ctr_0];
        const float xi_35 = _data_velocity_20_30_10[_stride_velocity_0 * ctr_0];
        const float xi_36 = _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0];
        const float xi_37 = _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0];
        const float xi_38 = _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0];
        const float xi_39 = _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0];
        const float xi_40 = _data_force_20_32_10[_stride_force_0 * ctr_0];
        const float xi_41 = _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0];
        const float xi_42 = _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0];
        const float xi_43 = _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0];
        const float xi_44 = _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0];
        const float xi_45 = _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0];
        const float xi_46 = _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0];
        const float xi_47 = _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0];
        const float xi_48 = _data_force_20_30_10[_stride_force_0 * ctr_0];
        const float xi_49 = _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0];
        const float xi_50 = _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0];
        const float xi_1 = xi_26;
        const float xi_2 = xi_27;
        const float xi_3 = xi_28;
        const float xi_4 = xi_29;
        const float xi_5 = xi_30;
        const float xi_6 = xi_31;
        const float xi_7 = xi_32;
        const float xi_8 = xi_33;
        const float xi_9 = xi_34;
        const float xi_10 = xi_35;
        const float xi_11 = xi_36;
        const float xi_12 = xi_37;
        const float xi_13 = xi_38;
        const float xi_14 = xi_39;
        const float xi_15 = xi_40;
        const float xi_16 = xi_41;
        const float xi_17 = xi_42;
        const float xi_18 = xi_43;
        const float xi_19 = xi_45;
        const float xi_20 = xi_46;
        const float xi_21 = xi_48;
        const float xi_22 = xi_47;
        const float xi_23 = xi_44;
        const float xi_24 = xi_49;
        const float xi_25 = xi_50;
        const float rho = xi_1 + xi_11 + xi_12 + xi_13 + xi_14 + xi_16 + xi_17 +
                          xi_18 + xi_19 + xi_2 + xi_20 + xi_22 + xi_23 + xi_24 +
                          xi_25 + xi_3 + xi_4 + xi_5 + xi_8;
        const float u1Pu2 = xi_7 + xi_9;
        const float u1Mu2 = -xi_7 + xi_9;
        _data_pdfs_20_30_10[_stride_pdfs_0 * ctr_0] =
            omega_shear * (rho * (xi_10 * xi_10) * -0.333333333333333f +
                           rho * (xi_7 * xi_7) * -0.333333333333333f +
                           rho * (xi_9 * xi_9) * -0.333333333333333f +
                           rho * 0.333333333333333f - xi_18) +
            xi_18 +
            (omega_shear * -0.5f + 1.0f) *
                (-xi_10 * xi_21 - xi_15 * xi_7 - xi_6 * xi_9);
        _data_pdfs_20_31_10[_stride_pdfs_0 * ctr_0] =
            omega_shear *
                (rho * (xi_7 * xi_7) * -0.166666666666667f +
                 rho * (xi_9 * xi_9) * -0.166666666666667f +
                 rho *
                     ((xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                   0.050000000000000003f) *
                      (xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                   0.050000000000000003f)) *
                     -0.166666666666667f +
                 rho * ((xi_9 * xi_9) * 0.333333333333333f +
                        xi_9 * 0.166666666666667f - 0.111111111111111f) +
                 rho * 0.166666666666667f - xi_25) +
            xi_25 +
            (omega_shear * -0.5f + 1.0f) *
                (xi_15 * xi_7 * -0.166666666666667f +
                 xi_21 *
                     (xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                  0.050000000000000003f) *
                     -0.166666666666667f +
                 xi_6 * (xi_9 * 2.0f + 1.0f) * 0.166666666666667f);
        _data_pdfs_20_32_10[_stride_pdfs_0 * ctr_0] =
            omega_shear *
                (rho * (xi_7 * xi_7) * -0.166666666666667f +
                 rho * (xi_9 * xi_9) * -0.166666666666667f +
                 rho *
                     ((xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                   0.050000000000000003f) *
                      (xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                   0.050000000000000003f)) *
                     -0.166666666666667f +
                 rho * ((xi_9 * xi_9) * 0.333333333333333f +
                        xi_9 * -0.166666666666667f - 0.111111111111111f) +
                 rho * 0.166666666666667f - xi_24) +
            xi_24 +
            (omega_shear * -0.5f + 1.0f) *
                (xi_15 * xi_7 * -0.166666666666667f +
                 xi_21 *
                     (xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                  0.050000000000000003f) *
                     -0.166666666666667f +
                 xi_6 * (xi_9 * 2.0f - 1.0f) * 0.166666666666667f);
        _data_pdfs_20_33_10[_stride_pdfs_0 * ctr_0] =
            omega_shear *
                (rho * (xi_10 * xi_10) * -0.166666666666667f +
                 rho * (xi_7 * xi_7) * -0.166666666666667f +
                 rho * (xi_9 * xi_9) * -0.166666666666667f +
                 rho * ((xi_10 * xi_10) * 0.333333333333333f +
                        xi_10 * -0.166666666666667f - 0.111111111111111f) +
                 rho * 0.166666666666667f - xi_17) +
            xi_17 +
            (omega_shear * -0.5f + 1.0f) *
                (xi_15 * xi_7 * -0.166666666666667f +
                 xi_21 * (xi_10 * 2.0f - 1.0f) * 0.166666666666667f +
                 xi_6 * xi_9 * -0.166666666666667f);
        _data_pdfs_20_34_10[_stride_pdfs_0 * ctr_0] =
            omega_shear *
                (rho * (xi_10 * xi_10) * -0.166666666666667f +
                 rho * (xi_7 * xi_7) * -0.166666666666667f +
                 rho * (xi_9 * xi_9) * -0.166666666666667f +
                 rho * ((xi_10 * xi_10) * 0.333333333333333f +
                        xi_10 * 0.166666666666667f - 0.111111111111111f) +
                 rho * 0.166666666666667f - xi_1) +
            xi_1 +
            (omega_shear * -0.5f + 1.0f) *
                (xi_15 * xi_7 * -0.166666666666667f +
                 xi_21 * (xi_10 * 2.0f + 1.0f) * 0.166666666666667f +
                 xi_6 * xi_9 * -0.166666666666667f);
        _data_pdfs_20_35_10[_stride_pdfs_0 * ctr_0] =
            omega_shear *
                (rho * (xi_10 * xi_10) * -0.166666666666667f +
                 rho * (xi_7 * xi_7) * -0.166666666666667f +
                 rho * (xi_9 * xi_9) * -0.166666666666667f +
                 rho * ((xi_7 * xi_7) * 0.333333333333333f +
                        xi_7 * 0.166666666666667f - 0.111111111111111f) +
                 rho * 0.166666666666667f - xi_3) +
            xi_3 +
            (omega_shear * -0.5f + 1.0f) *
                (xi_10 * xi_21 * -0.166666666666667f +
                 xi_15 * (xi_7 * 2.0f + 1.0f) * 0.166666666666667f +
                 xi_6 * xi_9 * -0.166666666666667f);
        _data_pdfs_20_36_10[_stride_pdfs_0 * ctr_0] =
            omega_shear *
                (rho * (xi_10 * xi_10) * -0.166666666666667f +
                 rho * (xi_7 * xi_7) * -0.166666666666667f +
                 rho * (xi_9 * xi_9) * -0.166666666666667f +
                 rho * ((xi_7 * xi_7) * 0.333333333333333f +
                        xi_7 * -0.166666666666667f - 0.111111111111111f) +
                 rho * 0.166666666666667f - xi_13) +
            xi_13 +
            (omega_shear * -0.5f + 1.0f) *
                (xi_10 * xi_21 * -0.166666666666667f +
                 xi_15 * (xi_7 * 2.0f - 1.0f) * 0.166666666666667f +
                 xi_6 * xi_9 * -0.166666666666667f);
        _data_pdfs_20_37_10[_stride_pdfs_0 * ctr_0] =
            omega_shear * (rho * (xi_7 * xi_7) * -0.0416666666666667f +
                           rho * (xi_9 * xi_9) * -0.0416666666666667f +
                           rho *
                               ((xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                             0.050000000000000003f) *
                                (xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                             0.050000000000000003f)) *
                               -0.0416666666666667f +
                           rho * (xi_10 * -0.0833333333333333f +
                                  (xi_7 * xi_7) * 0.0416666666666667f +
                                  xi_9 * 0.0833333333333333f +
                                  ((xi_10 - xi_9 +
                                    ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                        0.050000000000000003f) *
                                   (xi_10 - xi_9 +
                                    ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                        0.050000000000000003f)) *
                                      0.125f +
                                  ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                      -0.0041666666666666666f -
                                  0.0138888888888889f) +
                           rho * 0.0416666666666667f - xi_14) +
            xi_14 +
            (omega_shear * -0.5f + 1.0f) *
                (xi_15 * xi_7 * -0.0833333333333333f +
                 xi_21 *
                     (xi_10 * -2.0f + xi_9 * 3.0f +
                      ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                          -0.10000000000000001f +
                      1.0f) *
                     -0.0833333333333333f +
                 xi_6 *
                     (xi_10 * -3.0f + xi_9 * 2.0f +
                      ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                          -0.15000000000000002f +
                      1.0f) *
                     0.0833333333333333f);
        _data_pdfs_20_38_10[_stride_pdfs_0 * ctr_0] =
            omega_shear * (rho * (xi_7 * xi_7) * -0.0416666666666667f +
                           rho * (xi_9 * xi_9) * -0.0416666666666667f +
                           rho *
                               ((xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                             0.050000000000000003f) *
                                (xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                             0.050000000000000003f)) *
                               -0.0416666666666667f +
                           rho * (xi_10 * 0.0833333333333333f +
                                  (xi_7 * xi_7) * 0.0416666666666667f +
                                  xi_9 * 0.0833333333333333f +
                                  ((xi_10 + xi_9 +
                                    ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                        0.050000000000000003f) *
                                   (xi_10 + xi_9 +
                                    ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                        0.050000000000000003f)) *
                                      0.125f +
                                  ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                      0.0041666666666666666f -
                                  0.0138888888888889f) +
                           rho * 0.0416666666666667f - xi_12) +
            xi_12 +
            (omega_shear * -0.5f + 1.0f) *
                (xi_15 * xi_7 * -0.0833333333333333f +
                 xi_21 *
                     (xi_10 * 2.0f + xi_9 * 3.0f +
                      ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                          0.10000000000000001f +
                      1.0f) *
                     0.0833333333333333f +
                 xi_6 *
                     (xi_10 * 3.0f + xi_9 * 2.0f +
                      ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                          0.15000000000000002f +
                      1.0f) *
                     0.0833333333333333f);
        _data_pdfs_20_39_10[_stride_pdfs_0 * ctr_0] =
            omega_shear * (rho * (xi_7 * xi_7) * -0.0416666666666667f +
                           rho * (xi_9 * xi_9) * -0.0416666666666667f +
                           rho *
                               ((xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                             0.050000000000000003f) *
                                (xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                             0.050000000000000003f)) *
                               -0.0416666666666667f +
                           rho * (xi_10 * -0.0833333333333333f +
                                  (xi_7 * xi_7) * 0.0416666666666667f +
                                  xi_9 * -0.0833333333333333f +
                                  ((xi_10 + xi_9 +
                                    ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                        0.050000000000000003f) *
                                   (xi_10 + xi_9 +
                                    ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                        0.050000000000000003f)) *
                                      0.125f +
                                  ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                      -0.0041666666666666666f -
                                  0.0138888888888889f) +
                           rho * 0.0416666666666667f - xi_20) +
            xi_20 +
            (omega_shear * -0.5f + 1.0f) *
                (xi_15 * xi_7 * -0.0833333333333333f +
                 xi_21 *
                     (xi_10 * 2.0f + xi_9 * 3.0f +
                      ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.10000000000000001f -
                      1.0f) *
                     0.0833333333333333f +
                 xi_6 *
                     (xi_10 * 3.0f + xi_9 * 2.0f +
                      ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.15000000000000002f -
                      1.0f) *
                     0.0833333333333333f);
        _data_pdfs_20_310_10[_stride_pdfs_0 * ctr_0] =
            omega_shear * (rho * (xi_7 * xi_7) * -0.0416666666666667f +
                           rho * (xi_9 * xi_9) * -0.0416666666666667f +
                           rho *
                               ((xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                             0.050000000000000003f) *
                                (xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                             0.050000000000000003f)) *
                               -0.0416666666666667f +
                           rho * (xi_10 * 0.0833333333333333f +
                                  (xi_7 * xi_7) * 0.0416666666666667f +
                                  xi_9 * -0.0833333333333333f +
                                  ((xi_10 - xi_9 +
                                    ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                        0.050000000000000003f) *
                                   (xi_10 - xi_9 +
                                    ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                        0.050000000000000003f)) *
                                      0.125f +
                                  ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                      0.0041666666666666666f -
                                  0.0138888888888889f) +
                           rho * 0.0416666666666667f - xi_5) +
            xi_5 +
            (omega_shear * -0.5f + 1.0f) *
                (xi_15 * xi_7 * -0.0833333333333333f +
                 xi_21 *
                     (xi_10 * 2.0f + xi_9 * -3.0f +
                      ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.10000000000000001f +
                      1.0f) *
                     0.0833333333333333f +
                 xi_6 *
                     (xi_10 * 3.0f + xi_9 * -2.0f +
                      ((ctr_1 <= 0) ? (1.0f) : (0.0f)) * 0.15000000000000002f +
                      1.0f) *
                     -0.0833333333333333f);
        _data_pdfs_20_311_10[_stride_pdfs_0 * ctr_0] =
            omega_shear *
                (rho * (xi_7 * xi_7) * -0.0416666666666667f +
                 rho * (xi_9 * xi_9) * -0.0416666666666667f +
                 rho *
                     ((xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                   0.050000000000000003f) *
                      (xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                   0.050000000000000003f)) *
                     -0.0416666666666667f +
                 rho * ((u1Pu2 * u1Pu2) * 0.125f + u1Pu2 * 0.0833333333333333f +
                        ((xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                      0.050000000000000003f) *
                         (xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                      0.050000000000000003f)) *
                            0.0416666666666667f -
                        0.0138888888888889f) +
                 rho * 0.0416666666666667f - xi_23) +
            xi_23 +
            (omega_shear * -0.5f + 1.0f) *
                (xi_15 * (xi_7 * 2.0f + xi_9 * 3.0f + 1.0f) *
                     0.0833333333333333f +
                 xi_21 *
                     (xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                  0.050000000000000003f) *
                     -0.0833333333333333f +
                 xi_6 * (xi_7 * 3.0f + xi_9 * 2.0f + 1.0f) *
                     0.0833333333333333f);
        _data_pdfs_20_312_10[_stride_pdfs_0 * ctr_0] =
            omega_shear * (rho * (xi_7 * xi_7) * -0.0416666666666667f +
                           rho * (xi_9 * xi_9) * -0.0416666666666667f +
                           rho *
                               ((xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                             0.050000000000000003f) *
                                (xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                             0.050000000000000003f)) *
                               -0.0416666666666667f +
                           rho * ((u1Mu2 * u1Mu2) * 0.125f +
                                  u1Mu2 * -0.0833333333333333f +
                                  ((xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                                0.050000000000000003f) *
                                   (xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                                0.050000000000000003f)) *
                                      0.0416666666666667f -
                                  0.0138888888888889f) +
                           rho * 0.0416666666666667f - xi_16) +
            xi_16 +
            (omega_shear * -0.5f + 1.0f) *
                (xi_15 * (xi_7 * 2.0f + xi_9 * -3.0f + 1.0f) *
                     0.0833333333333333f +
                 xi_21 *
                     (xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                  0.050000000000000003f) *
                     -0.0833333333333333f +
                 xi_6 * (xi_7 * 3.0f + xi_9 * -2.0f + 1.0f) *
                     -0.0833333333333333f);
        _data_pdfs_20_313_10[_stride_pdfs_0 * ctr_0] =
            omega_shear * (rho * (xi_10 * xi_10) * -0.0416666666666667f +
                           rho * (xi_7 * xi_7) * -0.0416666666666667f +
                           rho * (xi_9 * xi_9) * -0.0416666666666667f +
                           rho * (xi_10 * -0.0833333333333333f +
                                  xi_7 * 0.0833333333333333f +
                                  (xi_9 * xi_9) * 0.0416666666666667f +
                                  ((xi_10 - xi_7) * (xi_10 - xi_7)) * 0.125f -
                                  0.0138888888888889f) +
                           rho * 0.0416666666666667f - xi_4) +
            xi_4 +
            (omega_shear * -0.5f + 1.0f) *
                (xi_15 * (xi_10 * -3.0f + xi_7 * 2.0f + 1.0f) *
                     0.0833333333333333f +
                 xi_21 * (xi_10 * -2.0f + xi_7 * 3.0f + 1.0f) *
                     -0.0833333333333333f +
                 xi_6 * xi_9 * -0.0833333333333333f);
        _data_pdfs_20_314_10[_stride_pdfs_0 * ctr_0] =
            omega_shear * (rho * (xi_10 * xi_10) * -0.0416666666666667f +
                           rho * (xi_7 * xi_7) * -0.0416666666666667f +
                           rho * (xi_9 * xi_9) * -0.0416666666666667f +
                           rho * (xi_10 * 0.0833333333333333f +
                                  xi_7 * 0.0833333333333333f +
                                  (xi_9 * xi_9) * 0.0416666666666667f +
                                  ((xi_10 + xi_7) * (xi_10 + xi_7)) * 0.125f -
                                  0.0138888888888889f) +
                           rho * 0.0416666666666667f - xi_22) +
            xi_22 +
            (omega_shear * -0.5f + 1.0f) *
                (xi_15 * (xi_10 * 3.0f + xi_7 * 2.0f + 1.0f) *
                     0.0833333333333333f +
                 xi_21 * (xi_10 * 2.0f + xi_7 * 3.0f + 1.0f) *
                     0.0833333333333333f +
                 xi_6 * xi_9 * -0.0833333333333333f);
        _data_pdfs_20_315_10[_stride_pdfs_0 * ctr_0] =
            omega_shear *
                (rho * (xi_7 * xi_7) * -0.0416666666666667f +
                 rho * (xi_9 * xi_9) * -0.0416666666666667f +
                 rho *
                     ((xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                   0.050000000000000003f) *
                      (xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                   0.050000000000000003f)) *
                     -0.0416666666666667f +
                 rho * ((u1Mu2 * u1Mu2) * 0.125f + u1Mu2 * 0.0833333333333333f +
                        ((xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                      0.050000000000000003f) *
                         (xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                      0.050000000000000003f)) *
                            0.0416666666666667f -
                        0.0138888888888889f) +
                 rho * 0.0416666666666667f - xi_19) +
            xi_19 +
            (omega_shear * -0.5f + 1.0f) *
                (xi_15 * (xi_7 * -2.0f + xi_9 * 3.0f + 1.0f) *
                     -0.0833333333333333f +
                 xi_21 *
                     (xi_10 + ((ctr_1 >= 63) ? (-1.0f) : (0.0f)) *
                                  0.050000000000000003f) *
                     -0.0833333333333333f +
                 xi_6 * (xi_7 * -3.0f + xi_9 * 2.0f + 1.0f) *
                     0.0833333333333333f);
        _data_pdfs_20_316_10[_stride_pdfs_0 * ctr_0] =
            omega_shear * (rho * (xi_7 * xi_7) * -0.0416666666666667f +
                           rho * (xi_9 * xi_9) * -0.0416666666666667f +
                           rho *
                               ((xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                             0.050000000000000003f) *
                                (xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                             0.050000000000000003f)) *
                               -0.0416666666666667f +
                           rho * ((u1Pu2 * u1Pu2) * 0.125f +
                                  u1Pu2 * -0.0833333333333333f +
                                  ((xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                                0.050000000000000003f) *
                                   (xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                                0.050000000000000003f)) *
                                      0.0416666666666667f -
                                  0.0138888888888889f) +
                           rho * 0.0416666666666667f - xi_8) +
            xi_8 +
            (omega_shear * -0.5f + 1.0f) *
                (xi_15 * (xi_7 * 2.0f + xi_9 * 3.0f - 1.0f) *
                     0.0833333333333333f +
                 xi_21 *
                     (xi_10 + ((ctr_1 <= 0) ? (1.0f) : (0.0f)) *
                                  0.050000000000000003f) *
                     -0.0833333333333333f +
                 xi_6 * (xi_7 * 3.0f + xi_9 * 2.0f - 1.0f) *
                     0.0833333333333333f);
        _data_pdfs_20_317_10[_stride_pdfs_0 * ctr_0] =
            omega_shear * (rho * (xi_10 * xi_10) * -0.0416666666666667f +
                           rho * (xi_7 * xi_7) * -0.0416666666666667f +
                           rho * (xi_9 * xi_9) * -0.0416666666666667f +
                           rho * (xi_10 * -0.0833333333333333f +
                                  xi_7 * -0.0833333333333333f +
                                  (xi_9 * xi_9) * 0.0416666666666667f +
                                  ((xi_10 + xi_7) * (xi_10 + xi_7)) * 0.125f -
                                  0.0138888888888889f) +
                           rho * 0.0416666666666667f - xi_11) +
            xi_11 +
            (omega_shear * -0.5f + 1.0f) *
                (xi_15 * (xi_10 * 3.0f + xi_7 * 2.0f - 1.0f) *
                     0.0833333333333333f +
                 xi_21 * (xi_10 * 2.0f + xi_7 * 3.0f - 1.0f) *
                     0.0833333333333333f +
                 xi_6 * xi_9 * -0.0833333333333333f);
        _data_pdfs_20_318_10[_stride_pdfs_0 * ctr_0] =
            omega_shear * (rho * (xi_10 * xi_10) * -0.0416666666666667f +
                           rho * (xi_7 * xi_7) * -0.0416666666666667f +
                           rho * (xi_9 * xi_9) * -0.0416666666666667f +
                           rho * (xi_10 * 0.0833333333333333f +
                                  xi_7 * -0.0833333333333333f +
                                  (xi_9 * xi_9) * 0.0416666666666667f +
                                  ((xi_10 - xi_7) * (xi_10 - xi_7)) * 0.125f -
                                  0.0138888888888889f) +
                           rho * 0.0416666666666667f - xi_2) +
            xi_2 +
            (omega_shear * -0.5f + 1.0f) *
                (xi_15 * (xi_10 * 3.0f + xi_7 * -2.0f + 1.0f) *
                     -0.0833333333333333f +
                 xi_21 * (xi_10 * 2.0f + xi_7 * -3.0f + 1.0f) *
                     0.0833333333333333f +
                 xi_6 * xi_9 * -0.0833333333333333f);
      }
    }
  }
}
} // namespace internal_ab1f3bc3368574afb482da84ccb58898

void CollideSweepSinglePrecisionLeesEdwards::run(IBlock *block) {
  auto pdfs = block->getData<field::GhostLayerField<float, 19>>(pdfsID);
  auto force = block->getData<field::GhostLayerField<float, 3>>(forceID);
  auto velocity = block->getData<field::GhostLayerField<float, 3>>(velocityID);

  auto &omega_shear = this->omega_shear_;
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(force->nrOfGhostLayers()));
  float *RESTRICT const _data_force = force->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(force->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(pdfs->nrOfGhostLayers()));
  float *RESTRICT _data_pdfs = pdfs->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(pdfs->layout(), field::fzyx);
  WALBERLA_ASSERT_GREATER_EQUAL(0, -int_c(velocity->nrOfGhostLayers()));
  float *RESTRICT const _data_velocity = velocity->dataAt(0, 0, 0, 0);
  WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx);
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
  const int64_t _stride_velocity_0 = int64_t(velocity->xStride());
  const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
  const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
  const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
  internal_ab1f3bc3368574afb482da84ccb58898::
      collidesweepsingleprecisionleesedwards_collidesweepsingleprecisionleesedwards(
          _data_force, _data_pdfs, _data_velocity, _size_force_0, _size_force_1,
          _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2,
          _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2,
          _stride_pdfs_3, _stride_velocity_0, _stride_velocity_1,
          _stride_velocity_2, _stride_velocity_3, omega_shear);
}

void CollideSweepSinglePrecisionLeesEdwards::runOnCellInterval(
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
  auto velocity = block->getData<field::GhostLayerField<float, 3>>(velocityID);

  auto &omega_shear = this->omega_shear_;
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
  WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin(), -int_c(velocity->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin(), -int_c(velocity->nrOfGhostLayers()));
  WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin(), -int_c(velocity->nrOfGhostLayers()));
  float *RESTRICT const _data_velocity =
      velocity->dataAt(ci.xMin(), ci.yMin(), ci.zMin(), 0);
  WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx);
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
  const int64_t _stride_velocity_0 = int64_t(velocity->xStride());
  const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
  const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
  const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
  internal_ab1f3bc3368574afb482da84ccb58898::
      collidesweepsingleprecisionleesedwards_collidesweepsingleprecisionleesedwards(
          _data_force, _data_pdfs, _data_velocity, _size_force_0, _size_force_1,
          _size_force_2, _stride_force_0, _stride_force_1, _stride_force_2,
          _stride_force_3, _stride_pdfs_0, _stride_pdfs_1, _stride_pdfs_2,
          _stride_pdfs_3, _stride_velocity_0, _stride_velocity_1,
          _stride_velocity_2, _stride_velocity_3, omega_shear);
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